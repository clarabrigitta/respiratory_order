# fit SEIRS model to data ----

# Selecting data for Rhinovirus Jan 2020 - Dec 2021 ----
data <- read_csv("inst/data/respiratory_scot_20250917.csv") %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>% 
  filter(!Pathogen %in% c("Influenza - Type A (not subtyped)",
                          "Influenza - Type A or B",
                          "Influenza - Type A(H1N1)pdm09" ,
                          "Influenza - Type A(H3)",
                          "Mycoplasma pneumoniae" ),
         WeekBeginning >= as.Date("2020-03-23"),
         WeekBeginning <= as.Date("2022-03-02"))

# function to run model, fitting omega
run_model <- function(I0, R0, omega, p_inf){
  N <- 5500000
  I0 <- I0
  E0 <- 0
  R0 <- R0
  S0 <- N - I0 - E0 - R0
  y0 <- c(S = S0, E = E0, I = I0, R = R0)
  
  out <- ode(y = y0, times = times, func = seirs, parms = list(N = N, sigma = sigma, gamma = gamma, omega = omega, p_inf = p_inf)) %>% 
    as.data.frame() %>% 
    mutate(date = as.Date(as.Date("2020-03-23"):as.Date("2022-03-02")), 
           week_beginning = floor_date(date, unit = "week", week_start = 1)) %>% 
    group_by(week_beginning) %>% 
    summarise(inc = sum(inc))
  
  return(out)
}

# bayesian model specification ----
likelihood <- function(param) {
  
  initial_I <- param[1]
  initial_R <- param[2]
  imm_period <- param[3]
  prob_inf <- param[4]
  detection_rate <- param[5]
  
  likelihood <- dpois(data$NumberCasesPerWeek,
                      run_model(I0 = initial_I, R0 = initial_R, omega = 1/imm_period, p_inf = prob_inf)$inc * detection_rate,
                      log = T)
  
  return(sum(likelihood)) 
}

prior <- createUniformPrior(lower = c(1, 1, 1, 0, 0), upper = c(2750000, 2750000, 720, 1, 1))

setup <- createBayesianSetup(likelihood = likelihood, prior = prior)

# run sampler ----
settings <- list(iterations = 50000, nrChains = 1, burnin = 25000, thin = 5)

set.seed(24)
out <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)

## retrieve results
getSample(out)
summary(out)
plot(out)
correlationPlot(out)

# run with samples ----

posterior <- getSample(out, start = 2)

traj <- lapply(1:nrow(posterior),
               function(r){
                 N <- 5500000
                 I0 <- as.numeric(posterior[r, 1])
                 E0 <- 0
                 R0 <- as.numeric(posterior[r, 2])
                 S0 <- N - I0 - E0 - R0
                 y0 <- c(S = S0, E = E0, I = I0, R = R0)
                 
                 ode(y = y0, times = times, func = seirs, parms = list(N = N, sigma = sigma, gamma = gamma, omega = posterior[r, 3], p_inf = posterior[r, 4])) %>% 
                   as.data.frame() %>% 
                   mutate(date = as.Date(as.Date("2020-03-23"):as.Date("2022-03-02")), 
                          week_beginning = floor_date(date, unit = "week", week_start = 1)) %>% 
                   group_by(week_beginning) %>% 
                   summarise(inc = sum(inc) * posterior[r, 5]) %>% 
                   select(inc)
               })

traj_all <- do.call(cbind, traj) %>% 
  as.data.frame() %>%
  bind_cols(data[, c(1, 4)]) %>% 
  pivot_longer(cols = c(1:length(traj)), names_to = "iter", values_to = "count")

traj_hdi <- test <- do.call(cbind, traj) %>% 
  as.data.frame() %>%
  t() %>% 
  hdi() %>%
  rbind(mean = colMeans(do.call(cbind, traj) %>% 
                          as.data.frame() %>%
                          t())) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(date = unique(data$WeekBeginning))


# plot trajectories ----
ggplot() + 
  geom_line(data = traj_all, aes(x = WeekBeginning, y = count, group = iter), alpha = 0.05, color = "darkgrey") +
  geom_line(data = data, aes(x = WeekBeginning, y = NumberCasesPerWeek), alpha = 0.7, color = "black") +
  ylim(0, 5000) +
  theme(legend.position = "none") +
  labs(x = "Months", y = "Count") + 
  theme_classic()

ggplot() +
  # geom_line(data = traj_hdi, aes(x = date, y = mean), color = "blue", size = 1) +
  # geom_ribbon(data = traj_hdi, aes(x = date, ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(data = combined, aes(x = date, y = mean_contacts*25, colour = part_age)) +
  geom_line(data = data, aes(x = WeekBeginning, y = NumberCasesPerWeek), alpha = 0.7, color = "red", size = 1) +
  scale_y_continuous(name = "Count", 
                     sec.axis = sec_axis(trans = ~./25, name = "Number of contacts")) +
  scale_color_viridis(discrete = T, option = "D") +
  labs(x = "Months", y = "Count") + 
  theme_classic() +
  facet_wrap(~Pathogen)
