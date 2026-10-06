source(here::here("R", "library_hpc.r"))

# define date for file saving
date <- format(Sys.Date(), "%d%m%Y")

# fixed inputs built by R/prepare_inputs.r
inputs <- readRDS(here("inst", "outdata", "hpc_inputs.rds"))
fortnight_matrix <- inputs$fortnight_matrix
scot_population  <- inputs$scot_population
scot_births      <- inputs$scot_births

source(here("R", "model_rcpp.r"))

# scenario to process: array task i -> i-th scenario (first one if run locally)
scenario_names <- sapply(combinations, function(x) x$scenario)
scenario <- unique(scenario_names)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]
combinations <- combinations[scenario_names == scenario]

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
  
  seirs_rcpp(y0 = y0,
             times = times,
             parms = list(sigma = sigma,
                          gamma = gamma,
                          omega = 1 / imm_days,
                          p_inf = p_inf,
                          bg    = bg,
                          mu_b  = daily_births,
                          mu_d  = hazard_death),
             contacts_prepped = contacts_prepped)
}

# aggregate daily incidence to 5 data bands by week
weekly_incidence <- function(ode_out) {
  daily_inc <- ode_out[, inc_cols, drop = FALSE]
  
  daily_agg <- cbind(daily_inc[, 1],
                     daily_inc[, 2] + daily_inc[, 3],
                     daily_inc[, 4] + daily_inc[, 5] + daily_inc[, 6],
                     daily_inc[, 7] + daily_inc[, 8],
                     daily_inc[, 9])
  
  rowsum(daily_agg, week_group)
}

# Rt calculation ----
S_cols <- seq(2, 34, by = 4)
C_by_day <- lapply(fortnight_lookup, function(f) fortnight_matrix[[f]]$matrix)

calculate_rt <- function(ode_out, p_inf, gamma) {
  S <- ode_out[, S_cols, drop = FALSE]
  N <- S + ode_out[, S_cols + 1] + ode_out[, S_cols + 2] + ode_out[, S_cols + 3]
  scale <- p_inf / gamma
  
  # K[i,j] = C[i,j] * S[i] / N[j] * p_inf / gamma
  rt_vec <- vapply(seq_len(nrow(ode_out)), function(d) {
    K <- scale * S[d, ] * t(t(C_by_day[[d]]) / N[d, ])
    max(Mod(eigen(K, symmetric = FALSE, only.values = TRUE)$values))
  }, numeric(1))
  
  data.frame(time = ode_out[, "time"], rt = rt_vec, date = dates_seq)
}

# run with posterior samples ----

results_traj <- list()
results_rt   <- list()

## load model fit output
results <- readRDS(file = here("inst", "outdata", date, paste0("parameters_", date, "_", scenario))) # change date if needed

for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  cat("Trajectory + Rt", n, ":", virus_name)
  
  posterior <- getSample(results[[n]], start = 3, thin = 1000)
  gamma <- 1 / combinations[[n]]$inf_period
  
  out <- lapply(seq_len(nrow(posterior)), function(r) {
    detection_rates <- posterior[r, idx$det]
    p_inf           <- posterior[r, idx$pinf]
    
    ode_out <- run_model_rcpp(p_sus_bands = posterior[r, idx$sus],
                              n           = n,
                              imm_days    = posterior[r, idx$imm],
                              imports     = 10^posterior[r, idx$imp],
                              p_inf       = p_inf)
    
    expected_out <- sweep(weekly_incidence(ode_out), 2, detection_rates, `*`)
    colnames(expected_out) <- c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")
    
    list(traj = expected_out,
         rt   = calculate_rt(ode_out, p_inf = p_inf, gamma = gamma))
  })
  
  results_traj[[virus_name]] <- lapply(out, `[[`, "traj")
  results_rt[[virus_name]]   <- lapply(out, `[[`, "rt")
  
  cat("  Done:", virus_name, "\n")
}

saveRDS(results_traj, file = here("inst", "outdata", date, paste0("traj_", date, "_", scenario))) # change date as needed
saveRDS(results_rt,   file = here("inst", "outdata", date, paste0("rt_", date, "_", scenario)))   # change date as needed
