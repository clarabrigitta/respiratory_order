source(here::here("R", "library_hpc.r"))

# define date for file saving
date <- format(Sys.Date(), "%d%m%Y")

# fixed inputs built by R/prepare_inputs.r
inputs <- readRDS(here("inst", "outdata", "hpc_inputs.rds"))
fortnight_matrix <- inputs$fortnight_matrix
scot_population  <- inputs$scot_population
scot_births      <- inputs$scot_births

source(here("R", "model_rcpp.r"))

# helpers ----
fit_start <- as.Date("2020-03-23")
fit_end   <- as.Date("2022-03-02")
times <- seq(0, as.integer(fit_end - fit_start), by = 1)
dates_seq  <- seq.Date(fit_start, fit_end, by = 1)
week_group <- as.character(floor_date(dates_seq, unit = "week", week_start = 1))

contacts_prepped <- prepare_contacts_cpp(fortnight_matrix, fortnight_lookup)

pathogen_map <- c(fluA = "Influenza A",
                  fluB = "Influenza B",
                  RSV  = "RSV",
                  hCOV = "Seasonal coronavirus",
                  AdV  = "Adenovirus",
                  RV   = "Rhinovirus",
                  hMPV = "HMPV",
                  PIV  = "Parainfluenza (Any Type)")

tots <- c(tot1, tot2, tot3, tot4, tot5, tot6, tot7, tot8, tot9)

# indexing for parameters (to keep better track)
idx <- list(det = 1:5, sus = 6:10, imm = 11, imp = 12, pinf = 13)
band_of_group <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)   # 9 model groups -> 5 data bands
u_idx <- c(idx$det, idx$sus, idx$imp, idx$pinf)
inc_cols <- 1 + 4 * 9 + seq_len(9)

# model runner ----

run_model_rcpp <- function(p_sus_bands, n, imm_days, imports = 1, p_inf) {
  sigma <- 1 / combinations[[n]]$inc_period
  gamma <- 1 / combinations[[n]]$inf_period
  bg    <- imports / sum(tots)
  
  R_init <- tots * (1 - p_sus_bands[band_of_group])
  S_init <- tots - R_init
  E_init <- bg * S_init / sigma
  I_init <- bg * S_init / gamma
  
  y0_mat <- rbind(S = S_init - E_init - I_init, E = E_init,
                  I = I_init, R = R_init)
  y0 <- setNames(as.vector(y0_mat),
                 paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4)))
  
  out <- seirs_rcpp(y0 = y0,
                    times = times,
                    parms = list(sigma = sigma,
                                 gamma = gamma,
                                 omega = 1 / imm_days,
                                 p_inf = p_inf,
                                 bg    = bg,
                                 mu_b  = daily_births,
                                 mu_d  = hazard_death),
                    contacts_prepped = contacts_prepped)
  
  daily_inc <- out[, inc_cols, drop = FALSE]
  
  daily_agg <- cbind(daily_inc[, 1],
                     daily_inc[, 2] + daily_inc[, 3],
                     daily_inc[, 4] + daily_inc[, 5] + daily_inc[, 6],
                     daily_inc[, 7] + daily_inc[, 8],
                     daily_inc[, 9])
  
  rowsum(daily_agg, week_group)
}

# run with posterior samples ----

results_traj <- list()

## load model fit output
results <- readRDS(file = here("inst", "outdata", paste0("parameters_", date))) # change date if needed

for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  pathogen_name <- pathogen_map[[virus_name]]
  cat("Trajectory", n, ":", virus_name)
  
  local({
    n_local <- n
    posterior <- getSample(results[[n_local]], start = 3, thin = 10000)
    
    traj <- lapply(seq_len(nrow(posterior)), function(r) {
      detection_rates <- posterior[r, idx$det]
      p_sus_bands     <- posterior[r, idx$sus]
      imm_duration    <- posterior[r, idx$imm]
      imports         <- 10^posterior[r, idx$imp]
      p_inf           <- posterior[r, idx$pinf]
      
      model_out <- run_model_rcpp(p_sus_bands, n_local, imm_duration, imports, p_inf)
      expected_out <- sweep(model_out, 2, detection_rates, `*`)
      colnames(expected_out) <- c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")
      expected_out
    })
    
    results_traj[[virus_name]] <<- traj
  })
  
  cat("  Done:", virus_name, "\n")
}

saveRDS(results_traj, file = here("inst", "outdata", paste0("traj_", date))) # change date as needed

# Rt calculation ----
## date data frame for assistance/reference
dates <- data.frame(date = seq(fit_start, fit_end, 1)) %>%
  mutate(time = 0:(n()-1),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date)) %>% 
  mutate(fortnight_n = as.integer(factor(fortnight)))

calculate_rt <- function(ode_out, p_inf, gamma, dates_df = dates) {
  t_vec <- ode_out[, "time"]
  
  rt_vec <- vapply(seq_along(t_vec), function(i) {
    t <- t_vec[i]
    
    S <- ode_out[i, seq(2, 34, by = 4)]
    E <- ode_out[i, seq(3, 35, by = 4)]
    I <- ode_out[i, seq(4, 36, by = 4)]
    R <- ode_out[i, seq(5, 37, by = 4)]
    N <- S + E + I + R
    
    C <- c_t(t)
    
    K <- outer(1:9, 1:9, function(i_id, j_id) {
      C[cbind(i_id, j_id)] * p_inf * S[i_id] / (N[j_id] * gamma)
    })
    
    Re(eigen(K, only.values = TRUE)$values[1])
  }, numeric(1))
  
  data.frame(time = t_vec, rt = rt_vec) %>%
    left_join(dates_df, by = "time")
}

results_rt <- list()

for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  pathogen_name <- pathogen_map[[virus_name]]
  cat("Rt trajectory", n, ":", virus_name)
  
  local({
    n_local <- n
    posterior <- getSample(results[[n_local]], start = 3, thin = 10000)
    
    rt_traj <- lapply(seq_len(nrow(posterior)), function(r) {
      p_sus_bands  <- posterior[r, idx$sus]
      imm_duration <- posterior[r, idx$imm]
      imports      <- 10^posterior[r, idx$imp]
      p_inf        <- posterior[r, idx$pinf]
      
      sigma <- 1 / combinations[[n_local]]$inc_period
      gamma <- 1 / combinations[[n_local]]$inf_period
      bg    <- imports / sum(tots)
      
      R_init <- tots * (1 - p_sus_bands[band_of_group])
      S_init <- tots - R_init
      E_init <- bg * S_init / sigma
      I_init <- bg * S_init / gamma
      
      y0_mat <- rbind(S = S_init - E_init - I_init, E = E_init,
                      I = I_init, R = R_init)
      y0 <- setNames(as.vector(y0_mat),
                     paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4)))
      
      ode_out <- seirs_rcpp(
        y0  = y0,
        times = times,
        parms = list(
          sigma = 1 / combinations[[n_local]]$inc_period,
          gamma = 1 / combinations[[n_local]]$inf_period,
          omega = 1 / imm_duration,
          p_inf = p_inf,
          bg = bg,
          mu_b = daily_births,
          mu_d = hazard_death
        ),
        contacts_prepped = contacts_prepped
      )
      
      calculate_rt(ode_out,
                   p_inf = p_inf,
                   gamma = 1 / combinations[[n_local]]$inf_period)
    })
    
    results_rt[[virus_name]] <<- rt_traj
  })
  
  cat("  Done:", virus_name, "\n")
}

saveRDS(results_rt, file = here("inst", "outdata", paste0("rt_", date))) # change date as needed
