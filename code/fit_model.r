# Selecting data for Rhinovirus Jan 2020 - Dec 2021 ----
data <- read_csv("data/respiratory_scot_20250917.csv") %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>% 
  filter(Pathogen == "Rhinovirus",
         WeekBeginning >= as.Date("2020-03-23"),
         WeekBeginning <= as.Date("2021-03-16"))

#  Helper functions ----
## Define contact function over time (contacts per person per day) based on CoMix
c_t <- function(t) {
  if (t < 74) {return(2.94)} else 
    if (t < 130) {return(3.71)} else
      if (t < 166) {return(5.17)} else    
        if (t < 219) {return(7.37)} else
          if (t < 256) {return(5.54)} else
            if (t < 273) {return(6.31)} else
              if (t < 287) {return(3.29)} else
                if (t < 352) {return(3.23)} else
                  return(5.27)}

## Convert day to month of year
day_to_month <- function(t) {
  d <- (floor(t) %% 365) + 1
  cum <- cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31))
  which(d <= cum)[1]
}

#  SEIR ODE (returns incidence) ----
seir <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    c <- c_t(t)         
    m <- month_foi[day_to_month(t)] # vector of monthly multiplier 
    beta_eff <- R_0 * gamma * c * m    
    lambda <- beta_eff * I / N
    
    dS <- -lambda * S
    dE <-  lambda * S - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I
    
    inc <- sigma * E
    
    list(c(dS, dE, dI, dR), inc = inc)
  })
}

#  Run SEIR model for Rhinovirus and aggregate weekly ----
## Fixed parameters
N   <- 1000000 # total population
sigma <- 1 / 3 # rate E to I
gamma <- 1 / 3 # rate I to R
R_0 <- 2

t_start <- min(data$WeekBeginning) # 2020-03-23 (first week of reporting after CoMix started)
t_end   <- max(data$WeekBeginning) + 6 # +6 so entire last week included (last week of reporting while CoMix running)
times   <- seq(0, as.integer(t_end - t_start), by = 1)

## Initial conditions
R0 <- 0
I0 <- 1
E0 <- 0
S0 <- N - I0 - E0 - R0
y0 <- c(S = S0, E = E0, I = I0, R = R0)

## Run model
run_model <- function(month_foi) {
  out <- ode(y = y0, times = times, func = seir, parms = list(N = N, sigma = sigma, gamma = gamma, R_0 = R_0, month_foi = month_foi)) %>% 
    as.data.frame() %>% 
    mutate(date = as.Date(as.Date("2020-03-23"):as.Date("2021-03-21")), 
           week_beginning = floor_date(date, unit = "week", week_start = 1)) %>% 
    group_by(week_beginning) %>% 
    summarise(inc = sum(inc)) # converting output to weekly to match surveillance data
  
  return(out)
}

# Bayesian model specification ----
likelihood <- function(theta) {
  log_m <- theta[1:12] # monthly multipliers on log scale
  m_raw <- exp(log_m)
  
  likelihood <- dpois(data$NumberCasesPerWeek,
                      run_model(month_foi = m_raw)$inc,
                      log = T)
  
  return(sum(likelihood)) 
}

lower <- c(rep(-3, 12))
upper <- c(rep( 3, 12))

prior <- createUniformPrior(lower = lower, upper = upper)

setup <- createBayesianSetup(likelihood = likelihood, prior = prior)

# Run sampler ----
settings <- list(iterations = 3000, nrChains = 1, burnin = 1000, thin = 5)

set.seed(24)
out <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)

## Retrieve results
getSample(out)
summary(out)
plot(out)

# Run with samples ----

posterior <- getSample(out)

traj <- lapply(1:nrow(posterior),
               function(r){
                 ode(y = y0, times = times, func = seir, parms = list(N = N, sigma = sigma, gamma = gamma, R_0 = R_0, month_foi = exp(posterior[r, ]))) %>% 
                   as.data.frame() %>% 
                   mutate(date = as.Date(as.Date("2020-03-23"):as.Date("2021-03-21")), 
                          week_beginning = floor_date(date, unit = "week", week_start = 1)) %>% 
                   group_by(week_beginning) %>% 
                   summarise(inc = sum(inc)) %>% 
                   select(inc)
               })

traj_all <- do.call(cbind, traj) %>% 
  as.data.frame() %>%
  setNames(as.character(1:405)) %>% 
  bind_cols(data[, 1:2]) %>% 
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
  

## Plot trajectories
ggplot() + 
  geom_line(data = traj, aes(x = WeekBeginning, y = count, group = iter), alpha = 0.05, color = "darkgrey") +
  geom_line(data = data, aes(x = WeekBeginning, y = NumberCasesPerWeek), alpha = 0.7, color = "black") +
  theme(legend.position = "none") +
  labs(x = "Months", y = "Count") + 
  theme_classic()

ggplot() +
  geom_line(data = traj_hdi, aes(x = date, y = mean), color = "blue") +
  geom_ribbon(data = traj_hdi, aes(x = date, ymin = lower, ymax = upper), alpha = 0.2)

## Calculate mean of monthly multipliers
posterior_hdi <- exp(hdi(posterior) %>% 
  rbind(mean = colMeans(posterior))) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(month = 1:12)
  
## Plot monthly multipliers
ggplot() +
  geom_line(data = posterior_hdi, aes(x = month, y = mean), color = "blue") +
  geom_ribbon(data = posterior_hdi, aes(x = month, ymin = lower, ymax = upper), alpha = 0.2)


              