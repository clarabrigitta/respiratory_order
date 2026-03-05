# fit SEIRS model to data ----

# Selecting data for RSV ----
data <- read_csv("inst/data/respiratory_scot_20250917.csv") %>%
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>%
  filter(!Pathogen %in% c("Influenza - Type A (not subtyped)",
                          "Influenza - Type A or B",
                          "Influenza - Type A(H1N1)pdm09" ,
                          "Influenza - Type A(H3)",
                          "Mycoplasma pneumoniae" ),
         WeekBeginning >= as.Date("2020-03-23"),
         WeekBeginning <= as.Date("2022-03-02")) %>%
  filter(Pathogen == "Respiratory syncytial virus") %>%
  head(101)

# function to run model; things to still fit: seeding timing
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

run_model <- function(mu, sigma){
  y0 <- c(S1 = tot1-5-floor(tot1*(1- imm_distribution(1, mu, sigma))), E1 = 0, I1 = 5, R1 = floor(tot1*(1- imm_distribution(1, mu, sigma))),
          S2 = tot2-5-floor(tot2*(1- imm_distribution(2, mu, sigma))), E2 = 0, I2 = 5, R2 = floor(tot2*(1- imm_distribution(2, mu, sigma))),
          S3 = tot3-5-floor(tot3*(1- imm_distribution(3, mu, sigma))), E3 = 0, I3 = 5, R3 = floor(tot3*(1- imm_distribution(3, mu, sigma))),
          S4 = tot4-5-floor(tot4*(1- imm_distribution(4, mu, sigma))), E4 = 0, I4 = 5, R4 = floor(tot4*(1- imm_distribution(4, mu, sigma))),
          S5 = tot5-5-floor(tot5*(1- imm_distribution(5, mu, sigma))), E5 = 0, I5 = 5, R5 = floor(tot5*(1- imm_distribution(5, mu, sigma))),
          S6 = tot6-5-floor(tot6*(1- imm_distribution(6, mu, sigma))), E6 = 0, I6 = 5, R6 = floor(tot6*(1- imm_distribution(6, mu, sigma))),
          S7 = tot7-5-floor(tot7*(1- imm_distribution(7, mu, sigma))), E7 = 0, I7 = 5, R7 = floor(tot7*(1- imm_distribution(7, mu, sigma))),
          S8 = tot8-5-floor(tot8*(1- imm_distribution(8, mu, sigma))), E8 = 0, I8 = 5, R8 = floor(tot8*(1- imm_distribution(8, mu, sigma))),
          S9 = tot9-5-floor(tot9*(1- imm_distribution(9, mu, sigma))), E9 = 0, I9 = 5, R9 = floor(tot9*(1- imm_distribution(9, mu, sigma))))

  out <- ode(y = y0, times = times, func = seirs,
             parms = list(sigma = 1 / combinations[[n]]$inc_period, gamma = 1 / combinations[[n]]$inf_period, omega = 1 / combinations[[n]]$imm_period, p_inf = combinations[[n]]$p_inf, mu_b = daily_births, mu_d = hazard_death)) %>%
    as.data.frame() %>%
    select(matches("inc|time")) %>%
    pivot_longer(cols = 2:10) %>%
    group_by(time) %>%
    summarise(inc = sum(value)) %>%
    mutate(date = as.Date(as.Date("2020-03-23"):as.Date("2022-02-26")),
           week_beginning = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(week_beginning) %>%
    summarise(inc = sum(inc))

  return(out)
}

run_model <- function(mu, sigma){
  # vectorise initial conditions: compute imm for all 9 groups in one call
  tots <- c(tot1, tot2, tot3, tot4, tot5, tot6, tot7, tot8, tot9)
  imm  <- imm_distribution(1:9, mu, sigma)
  R_init <- floor(tots * (1 - imm))
  y0_mat <- rbind(S = tots - 5 - R_init, E = rep(0, 9), I = rep(5, 9), R = R_init)
  y0 <- setNames(as.vector(y0_mat), paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4)))

  out_mat <- ode(y = y0, times = times, func = seirs,
                 parms = list(sigma = 1 / combinations[[n]]$inc_period, gamma = 1 / combinations[[n]]$inf_period,
                              omega = 1 / combinations[[n]]$imm_period, p_inf = combinations[[n]]$p_inf,
                              mu_b = daily_births, mu_d = hazard_death))

  # base R aggregation: avoids pivot_longer / group_by overhead on every likelihood call
  inc_cols  <- grep("^inc", colnames(out_mat))
  daily_inc <- rowSums(out_mat[, inc_cols, drop = FALSE])
  dates_seq <- seq.Date(as.Date("2020-03-23"), as.Date("2022-02-26"), by = 1)
  week_start <- floor_date(dates_seq, unit = "week", week_start = 1)
  weekly_inc <- tapply(daily_inc, week_start, sum)

  data.frame(week_beginning = as.Date(names(weekly_inc)),
             inc = as.numeric(weekly_inc))
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

likelihood <- function(param) {

  imm_mu <- param[1]
  imm_sigma <- param[2]
  detection_rate <- param[3]

  likelihood <- dpois(data$NumberCasesPerWeek,
                      run_model(mu = imm_mu, sigma = imm_sigma)$inc * detection_rate,
                      log = T)

  return(sum(likelihood))
}

prior <- createUniformPrior(lower = c(1, 1, 0), upper = c(10, 20, 1))

# setup <- createBayesianSetup(likelihood = likelihood, prior = prior)
# settings <- list(iterations = 5000, nrChains = 1, burnin = 2500, thin = 5)

# parallel = TRUE runs each chain on a separate core (DEzs requires >= 3 chains)
setup <- createBayesianSetup(likelihood = likelihood, prior = prior, parallel = TRUE)

# run sampler ----
settings <- list(iterations = 5000, nrChains = 1, burnin = 2000, thin = 5)

set.seed(24)
out <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)

## retrieve results
getSample(out)
summary(out)
plot(out)
# correlationPlot(out)