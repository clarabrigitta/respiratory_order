# run with posterior samples ----

results_traj <- list()

# Process outputs of model fits ----

## load combinations
source(here("R", "create_combinations.r"))

combinations <- create_combinations()

## load model fit output
results <- readRDS(file = here("inst", "outdata", "parameters_13072026")) # change date as needed

# generate model trajectories ----
## better if fit_model_rcpp.r run up to model runner
for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  pathogen_name <- pathogen_map[[virus_name]]
  cat("Trajectory", n, ":", virus_name)
  
  local({
    n_local <- n
    posterior <- getSample(results[[n_local]], start = 3, thin = 100)
    
    traj <- lapply(seq_len(nrow(posterior)), function(r) {
      detection_rates <- posterior[r, idx$det]
      p_sus_bands     <- posterior[r, idx$sus]
      imm_duration    <- posterior[r, idx$imm]
      imports         <- 10^posterior[r, idx$imp]
      
      model_out <- run_model_rcpp(p_sus_bands, n_local, imm_duration, imports)
      expected_out <- sweep(model_out, 2, detection_rates, `*`)
      colnames(expected_out) <- c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")
      expected_out
    })
    
    results_traj[[virus_name]] <<- traj
  })
  
  cat("  Done:", virus_name, "\n")
}

saveRDS(results_traj, file = here("inst", "outdata", "traj_12082026")) # change date as needed

# Rt calculation ----
## date data frame for assistance/reference
# dates <- data.frame(date = seq(as.Date("23-03-2020", format = "%d-%m-%Y"), as.Date("02-03-2022", format = "%d-%m-%Y"), 1)) %>%
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
    posterior <- getSample(results[[n_local]], start = 2, thin = 100)
    
    rt_traj <- lapply(seq_len(nrow(posterior)), function(r) {
      p_sus_bands  <- posterior[r, idx$sus]
      imm_duration <- posterior[r, idx$imm]
      imports      <- 10^posterior[r, idx$imp]
      
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
          p_inf = combinations[[n_local]]$p_inf,
          bg = bg,
          mu_b = daily_births,
          mu_d = hazard_death
        ),
        contacts_prepped = contacts_prepped
      )
      
      calculate_rt(ode_out,
                   p_inf = combinations[[n_local]]$p_inf,
                   gamma = 1 / combinations[[n_local]]$inf_period)
    })
    
    results_rt[[virus_name]] <<- rt_traj
  })
  
  cat("  Done:", virus_name, "\n")
}

saveRDS(results_rt, file = here("inst", "outdata", "rt_13072026")) # change date as needed
