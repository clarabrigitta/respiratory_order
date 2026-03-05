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

### load full dataset
data <- read_csv("inst/data/respiratory_scot_20250917.csv") %>%
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>%
  filter(!Pathogen %in% c("Influenza - Type A (not subtyped)",
                          "Influenza - Type A or B",
                          "Influenza - Type A(H1N1)pdm09",
                          "Influenza - Type A(H3)",
                          "Mycoplasma pneumoniae"),
         WeekBeginning >= as.Date("2020-03-23"),
         WeekBeginning <= as.Date("2022-03-02"))

### date constants used
dates_seq  <- seq.Date(as.Date("2020-03-23"), as.Date("2022-02-26"), by = 1)
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

  inc_cols  <- grep("^inc", colnames(out))
  daily_inc <- rowSums(out[, inc_cols, drop = FALSE])
  weekly_inc <- tapply(daily_inc, week_start, sum)

  data.frame(week_beginning = as.Date(names(weekly_inc)),
             inc = as.numeric(weekly_inc))
}

### bayesian specifications
settings <- list(iterations = 10000, burnin = 2000,  nrChains = 1)

### fit all combinations
results <- list()

for (n in seq_along(combinations)) {
  virus_name <- combinations[[n]]$name
  pathogen_name <- pathogen_map[[virus_name]]
  cat("Fitting combination", n, ":", virus_name)

  # filter for specific pathogen data
  subdata <- data %>%
    filter(Pathogen == pathogen_name) %>%
    head(101)

  local({
    n_local    <- n
    data_local <- subdata
    # specify priors based on which virus
    prior_local <- createUniformPrior(lower = combinations[[n_local]]$lb, upper = combinations[[n_local]]$ub)

    likelihood <- function(param) {
      sus_dist <- param[1]
      detection_rate <- param[2]
      imm_duration <- param[3]

      ll <- dpois(data_local$NumberCasesPerWeek,
                  run_model(p_sus = sus_dist, n = n_local, imm_days = imm_duration)$inc * detection_rate,
                  log = TRUE)
      return(sum(ll))
    }

    setup <- createBayesianSetup(likelihood = likelihood, prior = prior_local, parallel = FALSE, names = c("sus_dist", "detection_rate", "imm_duration"))

    set.seed(24)
    out <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)

    results[[virus_name]] <<- out
  })

  cat("  Done:", virus_name, "\n")
}

saveRDS(results, file = here("inst", "outdata", "parameters_sus_det_imm_05032026"))

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
    posterior <- getSample(results[[n_local]], start = 2)
    
    traj <- lapply(1:nrow(posterior),
                   function(r){
                     run_model(p_sus = posterior[r, 1], n = n_local, imm_days = posterior[r, 3])$inc * posterior[r, 2]
                   })
    
    results_traj[[virus_name]] <<- traj
  })
    
  cat("  Done:", virus_name, "\n")
}

saveRDS(results_traj, file = here("inst", "outdata", "traj_sus_det_imm_05032026"))

## process trajectories for plotting ----
### choose which virus trajectory and data
n = 8
virus_name <- combinations[[n]]$name
pathogen_name <- pathogen_map[[virus_name]]

traj <- results_traj[[virus_name]]
data <- all_data %>%
  filter(Pathogen == pathogen_name) %>%
  head(101)

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

ggplot() +
  geom_line(data = traj_hdi, aes(x = date, y = mean), color = "red", linewidth = 1) +
  geom_ribbon(data = traj_hdi, aes(x = date, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
  geom_line(data = data, aes(x = WeekBeginning, y = NumberCasesPerWeek), alpha = 0.7, color = "black", size = 1) +
  scale_x_date(date_breaks = "3 month") +
  labs(x = "Weeks", y = "Count", title = virus_name) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))