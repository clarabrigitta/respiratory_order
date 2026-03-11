# age-stratified SEIRS model ----

## helper data and functions ----
### date data frame for assistance/reference
dates <- data.frame(date = seq(as.Date("23-03-2020", format = "%d-%m-%Y"), as.Date("02-03-2022", format = "%d-%m-%Y"), 1)) %>% 
  mutate(time = 0:(n()-1),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date)) %>% 
  mutate(fortnight_n = as.integer(factor(fortnight)))

### define function to choose which fortnightly contact matrix to use 
fortnight_lookup <- dates$fortnight_n

c_t <- function(t) {
  fortnight_matrix[[fortnight_lookup[findInterval(t, dates$time)]]]$matrix
}

### create data frame defining model parameters
source(here("R", "create_combinations.R"))
combinations <- create_combinations()

### generate functional form of immunity given age (distribution) to use for how many people in R compartment for each age group (flat for now)
sus_distribution <- 0.9 # arbitrary value

### average population estimate between 2020-2022 for each age group
for (i in seq_len(nrow(scot_population))) {
  assign(paste0("tot", i), scot_population$average_n[i])
}

### define daily births and death rate
daily_births <- scot_births %>% select(births_daily) %>% pull() %>% mean() %>% floor() # load scot_births from explore_scotland.R
hazard_death <- daily_births/tot9 # for now keep death rate equal to birth rate, calculate per-capita death hazard

### define aging rates
a1 <- 1/(5*365) # 0-4
a2 <- 1/(6*365) # 5-10
a3 <- 1/(4*365) # 11-14
a4 <- 1/(10*365) # 15-24
a5 <- 1/(10*365) # 25-34
a6 <- 1/(10*365) # 35-44
a7 <- 1/(10*365) # 45-54
a8 <- 1/(10*365) # 55-64

# aging rate vector
a_vec <- c(a1, a2, a3, a4, a5, a6, a7, a8)

## model ODEs ----
seirs <- function(t, y, parms) {
  sigma <- parms$sigma
  gamma <- parms$gamma
  omega <- parms$omega
  p_inf <- parms$p_inf
  mu_b  <- parms$mu_b
  mu_d  <- parms$mu_d
  
  # extract compartments as vectors (9 for each class for the number of age-compartments)
  S <- y[seq(1, 33, by = 4)]
  E <- y[seq(2, 34, by = 4)]
  I <- y[seq(3, 35, by = 4)]
  R <- y[seq(4, 36, by = 4)]
  N <- S + E + I + R
  
  # force of infection (matrix multiply)
  contacts <- c_t(t)
  lambda <- as.vector(contacts %*% (I / N)) * p_inf
  
  # aging rates: groups 1-8 age out, group 9 loses to mortality
  a_out <- c(a_vec, mu_d)

  # influx from aging (births into group 1)
  S_in <- c(mu_b, a_vec * S[1:8])
  E_in <- c(0,    a_vec * E[1:8])
  I_in <- c(0,    a_vec * I[1:8])
  R_in <- c(0,    a_vec * R[1:8])
  
  # ODEs
  dS <- S_in - lambda * S + omega * R - a_out * S
  dE <- E_in + lambda * S - sigma * E  - a_out * E
  dI <- I_in + sigma * E  - gamma * I  - a_out * I
  dR <- R_in + gamma * I  - omega * R  - a_out * R
  
  # incidence
  inc <- sigma * E
  
  list(c(rbind(dS, dE, dI, dR)),
       inc1 = inc[1], inc2 = inc[2], inc3 = inc[3],
       inc4 = inc[4], inc5 = inc[5], inc6 = inc[6],
       inc7 = inc[7], inc8 = inc[8], inc9 = inc[9])
}

## run model ----

### starting compartment split for each age group
y0 <- c(S1 = tot1-5-floor(tot1*sus_distribution), E1 = 0, I1 = 5, R1 = floor(tot1*sus_distribution), 
        S2 = tot2-5-floor(tot2*sus_distribution), E2 = 0, I2 = 5, R2 = floor(tot2*sus_distribution),
        S3 = tot3-5-floor(tot3*sus_distribution), E3 = 0, I3 = 5, R3 = floor(tot3*sus_distribution),
        S4 = tot4-5-floor(tot4*sus_distribution), E4 = 0, I4 = 5, R4 = floor(tot4*sus_distribution),
        S5 = tot5-5-floor(tot5*sus_distribution), E5 = 0, I5 = 5, R5 = floor(tot5*sus_distribution),
        S6 = tot6-5-floor(tot6*sus_distribution), E6 = 0, I6 = 5, R6 = floor(tot6*sus_distribution),
        S7 = tot7-5-floor(tot7*sus_distribution), E7 = 0, I7 = 5, R7 = floor(tot7*sus_distribution),
        S8 = tot8-5-floor(tot8*sus_distribution), E8 = 0, I8 = 5, R8 = floor(tot8*sus_distribution),
        S9 = tot9-5-floor(tot9*sus_distribution), E9 = 0, I9 = 5, R9 = floor(tot9*sus_distribution))

### duration of model run
times <- seq(0, 709, by = 1)

### parallel run model
print(paste("start:", Sys.time()))
out <- mclapply(1:length(combinations),
                function(n) {
                  ode(y = y0, times = times, func = seirs, 
                      parms = list(sigma = 1 / combinations[[n]]$inc_period, gamma = 1 / combinations[[n]]$inf_period, omega = 1 / combinations[[n]]$imm_period, p_inf = combinations[[n]]$p_inf, mu_b = daily_births, mu_d = hazard_death))
                },
                mc.cores = 4)
print(paste("end:", Sys.time()))

### name list output
names(out) <- map_chr(combinations, 1)

## plot model output ----
### select for proper out/virus before plotting
fig <- out %>% 
  as.data.frame() %>% 
  rename(inc1 = inc1.E1, inc2 = inc2.E2, inc3 = inc3.E3, inc4 = inc4.E4, inc5 = inc5.E5, inc6 = inc6.E6, inc7 = inc7.E7, inc8 = inc8.E8, inc9 = inc9.E9) %>%
  pivot_longer(cols = 2:ncol(out), names_to = "compartment", values_to = "count") %>% 
  filter(grepl('inc', compartment)) %>% 
  left_join(dates, by = join_by(time)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = count, colour = compartment)) +
  scale_colour_viridis_d(option = "G", 
                         end = 0.8,
                         labels = c("inc1" = "0-4",
                                    "inc2" = "5-10",
                                    "inc3" = "11-14",
                                    "inc4" = "15-24",
                                    "inc5" = "25-34",
                                    "inc6" = "35-44",
                                    "inc7" = "45-44",
                                    "inc8" = "55-64",
                                    "inc9" = "65+")) +
  scale_x_date(date_breaks = "1 month") +
  labs(colour = "Age Group", y = "Count", x = "Time") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

### plot model output with observed data overlay
out_data <- out %>% 
  as.data.frame() %>% 
  rename(inc1 = inc1.E1, inc2 = inc2.E2, inc3 = inc3.E3, inc4 = inc4.E4, inc5 = inc5.E5, inc6 = inc6.E6, inc7 = inc7.E7, inc8 = inc8.E8, inc9 = inc9.E9) %>%
  pivot_longer(cols = 2:ncol(out), names_to = "compartment", values_to = "count") %>% 
  filter(grepl('inc', compartment)) %>% 
  left_join(dates, by = join_by(time))

obs_data <- read_csv("inst/data/respiratory_scot_20250917.csv") %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>% 
  filter(Pathogen == "Rhinovirus")

ggplot() +
  geom_line(data = out_data, aes(x = date, y = count, colour = compartment)) +
  geom_line(data = obs_data, aes(x = WeekBeginning, y = NumberCasesPerWeek), colour = "darkred") +
  scale_colour_viridis_d(option = "G", 
                         end = 0.8,
                         labels = c("inc1" = "0-4",
                                    "inc2" = "5-10",
                                    "inc3" = "11-14",
                                    "inc4" = "15-24",
                                    "inc5" = "25-34",
                                    "inc6" = "35-44",
                                    "inc7" = "45-44",
                                    "inc8" = "55-64",
                                    "inc9" = "65+")) +
  scale_x_date(    limits = as.Date(c("2020-01-01", "2023-01-01")), date_breaks = "1 month", date_labels = "%Y-%m") +
  labs(colour = "Age Group", y = "Count", x = "Time") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

### save plot
ggsave(filename = here("inst", "plots", "contact_fortnight_hMPV_130126.png"), plot = fig, width = 13, height = 9, dpi = 300)