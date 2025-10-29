## old contact values
# Helper function: Define contact function over time (contacts per person per day)
# c_t <- function(t) {
#   if (t < 74) {return(2.94)} else 
#     if (t < 130) {return(3.71)} else
#       if (t < 166) {return(5.17)} else    
#         if (t < 219) {return(7.37)} else
#           if (t < 256) {return(5.54)} else
#             if (t < 273) {return(6.31)} else
#               if (t < 287) {return(3.29)} else
#                 if (t < 352) {return(3.23)} else
#                   return(5.27)}

## function for contact values based on CoMix
c_t <- function(t) {
  if (t < 74) {return(2.60)} else
    if (t < 130) {return(3.56)} else
      if (t < 166) {return(5.25)} else
        if (t < 219) {return(6.72)} else
          if (t < 256) {return(5.16)} else
            if (t < 273) {return(6.02)} else
              if (t < 287) {return(3.04)} else
                if (t < 352) {return(2.91)} else
                  if (t < 360) {return(4.91)} else
                  return(5.06)}

## date dataframe for duration of CoMix
dates <- data.frame(date = seq(as.Date("23-03-2020", format = "%d-%m-%Y"), as.Date("02-03-2022", format = "%d-%m-%Y"), 1)) %>% 
  mutate(time = 1:n())

## SEIRS model
seirs <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      beta <- c_t(t)* p_inf # trying with 1/4 probability of infection
      
      dS <- - beta * S * (I/N) + omega * R
      dE <-  beta * S * (I/N) - sigma * E
      dI <-  sigma * E - gamma * I
      dR <-  gamma * I - omega * R
      
      inc <- sigma * E
      
      list(c(dS, dE, dI, dR), inc = inc)
    })
  }

### duration of periods
inc_period <- 4.98
inf_period <- 6.16
imm_period <- 358.9

### rates 
sigma <- 1 / inc_period # rate from E to I
gamma <- 1 / inf_period # rate from I to R
omega <- 1 / imm_period # rate from R to S

p_inf <- 0.0972

### starting population
N <- 5500000
I0 <- 1
E0 <- 0
R0 <- 1000
S0 <- N - I0 - E0 - R0
y0 <- c(S = S0, E = E0, I = I0, R = R0)

### duration of model run
times <- seq(0, 709, by = 1)

## run model
out <- ode(y = y0, times = times, func = seirs, parms = list(N = N, sigma = sigma, gamma = gamma, omega = omega, p_inf = p_inf))

## aggregate to weekly data
out <- as.data.frame(out) %>% 
  left_join(dates, join_by(time)) %>% 
  mutate(week_beginning = floor_date(date, unit = "week", week_start = 1)) %>% 
  group_by(week_beginning) %>% 
  summarise(inc = sum(inc))

## plot model output
out %>% 
  pivot_longer(cols = c(S, E, I, R, inc), names_to = "compartment", values_to = "count") %>% 
  mutate(compartment = factor(compartment, levels = c("S", "E", "I", "R", "inc"))) %>%
  # filter(compartment =="inc") %>%
  ggplot() +
  geom_line(aes(x = date, y = count, colour = compartment)) +
  # geom_line(data = pathogen_data %>% filter(!Pathogen %in% c("Rhinovirus", "Adenovirus")), aes(x = WeekBeginning, y = NumberCasesPerWeek, colour = Pathogen)) +
  scale_colour_viridis(discrete = TRUE, option = "D") +
  theme_bw(base_size = 10)
