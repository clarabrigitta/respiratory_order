# model fitting for all combinations ----

## helper data and functions ----
### pathogen name reference: name in combination data frame vs name in Scottish data
pathogen_map <- c(
  fluA = "Influenza - Type A (any subtype)",
  fluB = "Influenza - Type B",
  RSV = "Respiratory syncytial virus",
  hCOV = "Seasonal coronavirus (Non-SARS-CoV-2)",
  AdV = "Adenovirus",
  RV = "Rhinovirus",
  hMPV = "Human metapneumovirus",
  PIV = "Parainfluenza virus"
)
pathogen_map <- c(
  fluA = "Influenza A",
  fluB = "Influenza B",
  RSV = "RSV",
  hCOV = "Seasonal coronavirus",
  AdV = "Adenovirus",
  RV = "Rhinovirus",
  hMPV = "HMPV",
  PIV = "Parainfluenza (Any Type)"
)

### load full dataset
# data <- read_csv("inst/data/respiratory_scot_20250917.csv") %>%
#   mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>%
#   filter(!Pathogen %in% c("Influenza - Type A (not subtyped)",
#                           "Influenza - Type A or B",
#                           "Influenza - Type A(H1N1)pdm09",
#                           "Influenza - Type A(H3)",
#                           "Mycoplasma pneumoniae"),
#          WeekBeginning >= as.Date("2020-03-23"),
#          WeekBeginning <= as.Date("2022-03-02"))

data <- read_csv("inst/data/cases_all_respiratory_pathogens_by_agegroup_sex_20251008.csv") %>% 
  filter(Sex == "Total", !AgeGroup %in% c("Total", "Unknown"), !Pathogen %in% c("Influenza (All)", "COVID-19", "Mycoplasma pneumoniae")) %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>% 
  filter(WeekBeginning >= as.Date("2020-03-23"),
         WeekBeginning <= as.Date("2022-03-02")) %>% 
  pivot_wider(names_from = AgeGroup, values_from = NumberCasesPerWeek) %>% 
  mutate(`0 to 4` = `<1` + `1 to 4`,
         `65+` = `65 to 74` + `75+`) %>% 
  select(-c(`<1`, `1 to 4`, `65 to 74`, `75+`)) %>% 
  arrange(WeekBeginning) %>% 
  select(Pathogen, "0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")

### date constants used
dates_seq  <- seq.Date(as.Date("2020-03-23"), as.Date("2022-03-02"), by = 1)
week_start <- floor_date(dates_seq, unit = "week", week_start = 1)

### population totals by age group
tots <- c(tot1, tot2, tot3, tot4, tot5, tot6, tot7, tot8, tot9)

## run model fitting ----
run_model <- function(p_sus, n, imm_days) {
  R_init <- floor(tots * (1 - p_sus))
  y0_mat <- rbind(S = tots - 5 - R_init, E = rep(0, 9), I = rep(5, 9), R = R_init)
  y0 <- setNames(as.vector(y0_mat), paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4)))

  out <- ode(y = y0, times = times, func = seirs,
             parms = list(sigma = 1 / combinations[[n]]$inc_period, gamma = 1 / combinations[[n]]$inf_period,
                          omega = 1 / imm_days, p_inf = combinations[[n]]$p_inf,
                          mu_b = daily_births, mu_d = hazard_death))

  # inc_cols  <- grep("^inc", colnames(out))
  # daily_inc <- rowSums(out[, inc_cols, drop = FALSE])
  # weekly_inc <- tapply(daily_inc, week_start, sum)
  inc_cols  <- grep("^inc", colnames(out))
  daily_inc <- out[, inc_cols, drop = FALSE]
  daily_agg <- cbind(daily_inc[, 1],
                         daily_inc[, 2] + daily_inc[, 3],
                         daily_inc[, 4] + daily_inc[, 5] + daily_inc[, 6],
                         daily_inc[, 7] + daily_inc[, 8],
                         daily_inc[, 9])
  weekly_inc <- rowsum(daily_agg, week_start)
  # data.frame(week_beginning = as.Date(names(weekly_inc)),
  #            inc = as.numeric(weekly_inc))
  return(weekly_inc)
}

### bayesian specifications
settings <- list(iterations = 100000, burnin = 50000,  nrChains = 1)

### fit all combinations
results <- list()

for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  pathogen_name <- pathogen_map[[virus_name]]
  cat("Fitting combination", n, ":", virus_name)

  # filter for specific pathogen data
  subdata <- data %>%
    filter(Pathogen == pathogen_name)
  
  local({
    n_local    <- n
    data_local <- subdata
    # specify priors based on which virus
    prior_local <- createUniformPrior(lower = combinations[[n_local]]$lb, upper = combinations[[n_local]]$ub)

    likelihood <- function(param) {
      detection_rates <- param[1:5]
      sus_dist <- param[6]
      imm_duration <- param[7]
      
      model_out <- run_model(p_sus = sus_dist, n = n_local, imm_days = imm_duration)
      expected_out <- sweep(model_out, 2, detection_rates, `*`)

      # ll <- dpois(data_local$NumberCasesPerWeek,
      #             run_model(p_sus = sus_dist, n = n_local, imm_days = imm_duration)$inc * detection_rate,
      #             log = TRUE)
      
      ll <- dpois(as.vector(as.matrix(data_local[, -1])), as.vector(expected_out), log = TRUE)
      return(sum(ll, na.rm = TRUE))
      
      return(sum(ll))
    }

    setup <- createBayesianSetup(likelihood = likelihood, prior = prior_local, parallel = FALSE, 
                                 names = c("detection_0to4", "detection_5to14", "detection_15to44", "detection_45to64", "detection_65plus", "sus_dist", "imm_duration"))

    set.seed(24)
    out <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)

    results[[virus_name]] <<- out
  })

  cat("  Done:", virus_name, "\n")
}

saveRDS(results, file = here("inst", "outdata", "parameters_sus_det_imm_16032026"))

### manual diagnostic checks
summary(results[["RV"]])
getSample(results[["PIV"]], start = 2)
plot(results[["RV"]], start = 2)
correlationPlot(results[["fluA"]])

## run with posterior samples ----
results_traj <- list()

for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  pathogen_name <- pathogen_map[[virus_name]]
  cat("Trajectory", n, ":", virus_name)
  
  local({
    n_local   <- n
    posterior <- getSample(results[[n_local]], start = 10)
    
    traj <- lapply(1:nrow(posterior),
                   function(r){
                     detection_rates <- posterior[r, 1:5]
                     model_out <- run_model(p_sus = posterior[r, 6], n = n_local, imm_days = posterior[r, 7])
                     expected_out <- sweep(model_out, 2, detection_rates, `*`)
                     colnames(expected_out) <- c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")
                     return(expected_out)
                   })
    
    results_traj[[virus_name]] <<- traj
  })
    
  cat("  Done:", virus_name, "\n")
}


saveRDS(results_traj, file = here("inst", "outdata", "traj_sus_det_imm_16032026"))

## process trajectories for plotting ----
### choose which virus trajectory and data
n = 6
virus_name <- combinations[[n]]$name
pathogen_name <- pathogen_map[[virus_name]]

traj <- results_traj[[virus_name]]
data <- all_data %>%
  filter(Pathogen == pathogen_name) 

traj_all <- do.call(cbind, traj) %>%
  as.data.frame() %>%
  bind_cols(data[, c(1, 2)]) %>%
  pivot_longer(cols = c(1:length(traj)), names_to = "iter", values_to = "count")

traj_hdi <- do.call(cbind, traj) %>%
  as.data.frame() %>%
  t() %>%
  hdi() %>%
  rbind(mean = colMeans(do.call(cbind, traj) %>%
                          as.data.frame() %>%
                          t())) %>%
  t() %>%
  as.data.frame() %>%
  mutate(date = unique(data$WeekBeginning))

## plot trajectories ----
ggplot() +
  geom_line(data = traj_all, aes(x = WeekBeginning, y = count, group = iter), alpha = 0.05, color = "darkgrey") +
  geom_line(data = data, aes(x = WeekBeginning, y = NumberCasesPerWeek), alpha = 0.7, color = "black") +
  # ylim(0, 5000) +
  theme(legend.position = "none") +
  labs(x = "Weeks", y = "Count") +
  theme_classic()

scale_factor <- max(data$NumberCasesPerWeek, na.rm = TRUE) / max(combined$mean_contacts, na.rm = TRUE)

ggplot() +
  geom_line(data = traj_hdi, aes(x = date, y = mean), color =  "red", linewidth = 1) +
  geom_ribbon(data = traj_hdi, aes(x = date, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  geom_line(data = data, aes(x = WeekBeginning, y = NumberCasesPerWeek), alpha = 0.9, color = "black", size = 1) +
  geom_line(data = combined, aes(x = date, y = mean_contacts * scale_factor), colour = "blue", lty = 2, size = 1) +
  scale_x_date(date_breaks = "3 month") +
  scale_y_continuous(name = "Count",
                     sec.axis = sec_axis(~ . / scale_factor, name = "Mean number of contacts", breaks = seq(0, ceiling(max(combined$mean_contacts)), by = 2))) +
  labs(x = "Weeks", title = virus_name) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        axis.title.y.right = element_text(colour = "blue"),
        axis.text.y.right = element_text(colour = "blue"),
        axis.line.y.right = element_line(colour = "blue"),
        axis.ticks.y.right = element_line(colour = "blue"))

## Rt from ODE output ----
calculate_rt <- function(ode_out, p_inf, gamma, dates_df = dates) {
  times <- ode_out[, "time"]
  
  rt_vec <- vapply(seq_along(times), function(i) {
    t <- times[i]
    
    # Extract S and N for each age group
    S <- ode_out[i, seq(2, 34, by = 4)]
    E <- ode_out[i, seq(3, 35, by = 4)]
    I <- ode_out[i, seq(4, 36, by = 4)]
    R <- ode_out[i, seq(5, 37, by = 4)]
    N <- S + E + I + R
    
    # Contact matrix at time t
    C <- c_t(t)
    
    # Next-generation matrix: K_ij = C_ij * p_inf * S_i / (N_j * gamma)
    K <- outer(1:9, 1:9, function(i_id, j_id) {
      C[cbind(i_id, j_id)] * p_inf * S[i_id] / (N[j_id] * gamma)
    })
    
    # Rt = spectral radius of K
    Re(eigen(K, only.values = TRUE)$values[1])
  }, numeric(1))
  
  data.frame(time = times, rt = rt_vec) %>%
    left_join(dates_df, by = "time")
}

### run for each virus
results_rt <- list()

for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  pathogen_name <- pathogen_map[[virus_name]]
  cat("Trajectory", n, ":", virus_name)
  
  local({
    n_local   <- n
    posterior <- getSample(results[[n_local]], start = 2, thin = 100)
    
    rt_traj <- lapply(1:nrow(posterior), function(r) {
      p_sus      <- posterior[r, 1]
      imm_days   <- posterior[r, 3]
      
      R_init <- floor(tots * (1 - p_sus))
      y0_mat <- rbind(S = tots - 5 - R_init, E = rep(0, 9), I = rep(5, 9), R = R_init)
      y0 <- setNames(as.vector(y0_mat), paste0(rep(c("S","E","I","R"), 9), rep(1:9, each=4)))
      
      ode_out <- ode(y = y0, times = times, func = seirs,
                     parms = list(sigma = 1/combinations[[n_local]]$inc_period,
                                  gamma = 1/combinations[[n_local]]$inf_period,
                                  omega = 1/imm_days,
                                  p_inf = combinations[[n_local]]$p_inf,
                                  mu_b  = daily_births,
                                  mu_d  = hazard_death))
      
      calculate_rt(ode_out,
                   p_inf  = combinations[[n_local]]$p_inf,
                   gamma  = 1/combinations[[n_local]]$inf_period)
    })
    
    results_rt[[virus_name]] <<- rt_traj
  })
  
  cat("  Done:", virus_name, "\n")
}

saveRDS(results_rt, file = here("inst", "outdata", "rt_sus_det_imm_09032026"))
