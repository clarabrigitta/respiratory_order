# fit SEIRS model to data ----

# function to run model, fitting omega
run_model <- function(omega){
  out <- ode(y = y0, times = times, func = seirs, parms = list(N = N, sigma = sigma, gamma = gamma, omega = omega)) %>% 
    as.data.frame() %>% 
    mutate(date = as.Date(as.Date("2020-03-23"):as.Date("2022-03-02")), 
           week_beginning = floor_date(date, unit = "week", week_start = 1)) %>% 
    group_by(week_beginning) %>% 
    summarise(inc = sum(inc))
  
  return(out)
}

# bayesian model specification ----
likelihood <- function(param) {
  
  imm_period <- param[1]
  detection_rate <- param[2]
  
  likelihood <- dpois(data$NumberCasesPerWeek,
                      run_model(omega = 1/imm_period)$inc * detection_rate,
                      log = T)
  
  return(sum(likelihood)) 
}

prior <- createUniformPrior(lower = c(1, 0), upper = c(730, 1))

setup <- createBayesianSetup(likelihood = likelihood, prior = prior)

# run sampler ----
settings <- list(iterations = 3000, nrChains = 1, burnin = 1000, thin = 5)

set.seed(24)
out <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)

## retrieve results
getSample(out)
summary(out)
plot(out)

# run with samples ----

posterior <- getSample(out)

traj <- lapply(1:nrow(posterior),
               function(r){
                 ode(y = y0, times = times, func = seirs, parms = list(N = N, sigma = sigma, gamma = gamma, omega = posterior[r, 1])) %>% 
                   as.data.frame() %>% 
                   mutate(date = as.Date(as.Date("2020-03-23"):as.Date("2022-03-02")), 
                          week_beginning = floor_date(date, unit = "week", week_start = 1)) %>% 
                   group_by(week_beginning) %>% 
                   summarise(inc = sum(inc) * posterior[r, 2]) %>% 
                   select(inc)
               })

traj_all <- do.call(cbind, traj) %>% 
  as.data.frame() %>%
  setNames(as.character(1:405)) %>% 
  bind_cols(data[, c(1, 4)]) %>% 
  pivot_longer(cols = c(1:length(traj)), names_to = "iter", values_to = "count")

traj_hdi <- test <- do.call(cbind, traj) %>% 
  as.data.frame() %>%
  setNames(as.character(1:405)) %>% 
  t() %>% 
  hdi() %>%
  rbind(mean = colMeans(do.call(cbind, traj) %>% 
                          as.data.frame() %>%
                          setNames(as.character(1:405)) %>% 
                          t())) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(date = unique(data$WeekBeginning))


# plot trajectories ----
ggplot() + 
  geom_line(data = traj_all, aes(x = WeekBeginning, y = count, group = iter), alpha = 0.05, color = "darkgrey") +
  geom_line(data = data, aes(x = WeekBeginning, y = NumberCasesPerWeek), alpha = 0.7, color = "black") +
  theme(legend.position = "none") +
  ylim(0, 5000) +
  labs(x = "Months", y = "Count") + 
  theme_classic()

ggplot() +
  geom_line(data = traj_hdi, aes(x = date, y = mean), color = "blue") +
  geom_ribbon(data = traj_hdi, aes(x = date, ymin = lower, ymax = upper), alpha = 0.2)