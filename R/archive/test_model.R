### define function to choose which fortnightly contact matrix to use (choose the contact matrix based on a selection from fortnight_matrix)
c_t <- function(t){
  #filter fortnight_n value/row based on t, t defined in ode but is continuous and not discrete, using findInterval to map to discrete
  fortnight_value <- filter(dates, time == findInterval(t, dates$time)) %>%
    select(fortnight_n) %>%
    as.integer()

  #return correct fortnight_matrix contact matrix based on fortnight_n value
  return(fortnight_matrix[[fortnight_value]]$matrix)
}

## model ODEs ----
seirs <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      S1 <- y[1]
      E1 <- y[2]
      I1 <- y[3]
      R1 <- y[4]
      S2 <- y[5]
      E2 <- y[6]
      I2 <- y[7]
      R2 <- y[8]
      S3 <- y[9]
      E3 <- y[10]
      I3 <- y[11]
      R3 <- y[12]
      S4 <- y[13]
      E4 <- y[14]
      I4 <- y[15]
      R4 <- y[16]
      S5 <- y[17]
      E5 <- y[18]
      I5 <- y[19]
      R5 <- y[20]
      S6 <- y[21]
      E6 <- y[22]
      I6 <- y[23]
      R6 <- y[24]
      S7 <- y[25]
      E7 <- y[26]
      I7 <- y[27]
      R7 <- y[28]
      S8 <- y[29]
      E8 <- y[30]
      I8 <- y[31]
      R8 <- y[32]
      S9 <- y[33]
      E9 <- y[34]
      I9 <- y[35]
      R9 <- y[36]

      N1 <- S1 + E1 + I1 + R1
      N2 <- S2 + E2 + I2 + R2
      N3 <- S3 + E3 + I3 + R3
      N4 <- S4 + E4 + I4 + R4
      N5 <- S5 + E5 + I5 + R5
      N6 <- S6 + E6 + I6 + R6
      N7 <- S7 + E7 + I7 + R7
      N8 <- S8 + E8 + I8 + R8
      N9 <- S9 + E9 + I9 + R9

      contacts <- c_t(t)

      # force of infection calculations for 9 age groups
      lambda1 <- (contacts[1,1] * p_inf) * (I1/N1) + (contacts[1,2] * p_inf) * (I2/N2) + (contacts[1,3] * p_inf) * (I3/N3) + (contacts[1,4] * p_inf) * (I4/N4) + (contacts[1,5] * p_inf) * (I5/N5) + (contacts[1,6] * p_inf) * (I6/N6) + (contacts[1,7] * p_inf) * (I7/N7) + (contacts[1,8] * p_inf) * (I8/N8) + (contacts[1,9] * p_inf) * (I9/N9)
      lambda2 <- (contacts[2,1] * p_inf) * (I1/N1) + (contacts[2,2] * p_inf) * (I2/N2) + (contacts[2,3] * p_inf) * (I3/N3) + (contacts[2,4] * p_inf) * (I4/N4) + (contacts[2,5] * p_inf) * (I5/N5) + (contacts[2,6] * p_inf) * (I6/N6) + (contacts[2,7] * p_inf) * (I7/N7) + (contacts[2,8] * p_inf) * (I8/N8) + (contacts[2,9] * p_inf) * (I9/N9)
      lambda3 <- (contacts[3,1] * p_inf) * (I1/N1) + (contacts[3,2] * p_inf) * (I2/N2) + (contacts[3,3] * p_inf) * (I3/N3) + (contacts[3,4] * p_inf) * (I4/N4) + (contacts[3,5] * p_inf) * (I5/N5) + (contacts[3,6] * p_inf) * (I6/N6) + (contacts[3,7] * p_inf) * (I7/N7) + (contacts[3,8] * p_inf) * (I8/N8) + (contacts[3,9] * p_inf) * (I9/N9)
      lambda4 <- (contacts[4,1] * p_inf) * (I1/N1) + (contacts[4,2] * p_inf) * (I2/N2) + (contacts[4,3] * p_inf) * (I3/N3) + (contacts[4,4] * p_inf) * (I4/N4) + (contacts[4,5] * p_inf) * (I5/N5) + (contacts[4,6] * p_inf) * (I6/N6) + (contacts[4,7] * p_inf) * (I7/N7) + (contacts[4,8] * p_inf) * (I8/N8) + (contacts[4,9] * p_inf) * (I9/N9)
      lambda5 <- (contacts[5,1] * p_inf) * (I1/N1) + (contacts[5,2] * p_inf) * (I2/N2) + (contacts[5,3] * p_inf) * (I3/N3) + (contacts[5,4] * p_inf) * (I4/N4) + (contacts[5,5] * p_inf) * (I5/N5) + (contacts[5,6] * p_inf) * (I6/N6) + (contacts[5,7] * p_inf) * (I7/N7) + (contacts[5,8] * p_inf) * (I8/N8) + (contacts[5,9] * p_inf) * (I9/N9)
      lambda6 <- (contacts[6,1] * p_inf) * (I1/N1) + (contacts[6,2] * p_inf) * (I2/N2) + (contacts[6,3] * p_inf) * (I3/N3) + (contacts[6,4] * p_inf) * (I4/N4) + (contacts[6,5] * p_inf) * (I5/N5) + (contacts[6,6] * p_inf) * (I6/N6) + (contacts[6,7] * p_inf) * (I7/N7) + (contacts[6,8] * p_inf) * (I8/N8) + (contacts[6,9] * p_inf) * (I9/N9)
      lambda7 <- (contacts[7,1] * p_inf) * (I1/N1) + (contacts[7,2] * p_inf) * (I2/N2) + (contacts[7,3] * p_inf) * (I3/N3) + (contacts[7,4] * p_inf) * (I4/N4) + (contacts[7,5] * p_inf) * (I5/N5) + (contacts[7,6] * p_inf) * (I6/N6) + (contacts[7,7] * p_inf) * (I7/N7) + (contacts[7,8] * p_inf) * (I8/N8) + (contacts[7,9] * p_inf) * (I9/N9)
      lambda8 <- (contacts[8,1] * p_inf) * (I1/N1) + (contacts[8,2] * p_inf) * (I2/N2) + (contacts[8,3] * p_inf) * (I3/N3) + (contacts[8,4] * p_inf) * (I4/N4) + (contacts[8,5] * p_inf) * (I5/N5) + (contacts[8,6] * p_inf) * (I6/N6) + (contacts[8,7] * p_inf) * (I7/N7) + (contacts[8,8] * p_inf) * (I8/N8) + (contacts[8,9] * p_inf) * (I9/N9)
      lambda9 <- (contacts[9,1] * p_inf) * (I1/N1) + (contacts[9,2] * p_inf) * (I2/N2) + (contacts[9,3] * p_inf) * (I3/N3) + (contacts[9,4] * p_inf) * (I4/N4) + (contacts[9,5] * p_inf) * (I5/N5) + (contacts[9,6] * p_inf) * (I6/N6) + (contacts[9,7] * p_inf) * (I7/N7) + (contacts[9,8] * p_inf) * (I8/N8) + (contacts[9,9] * p_inf) * (I9/N9)

      # ODEs for 9 age groups
      dS1 <- mu_b - lambda1 * S1 + omega * R1 - a1 * S1
      dE1 <-  lambda1 * S1 - sigma * E1 - a1 * E1
      dI1 <-  sigma * E1 - gamma * I1 - a1 * I1
      dR1 <-  gamma * I1 - omega * R1 - a1 * R1

      dS2 <- - lambda2 * S2 + omega * R2 + a1 * S1 - a2 * S2
      dE2 <-  lambda2 * S2 - sigma * E2 + a1 * E1 - a2 * E2
      dI2 <-  sigma * E2 - gamma * I2 + a1 * I1 - a2 * I2
      dR2 <-  gamma * I2 - omega * R2 + a1 * R1 - a2 * R2

      dS3 <- - lambda3 * S3 + omega * R3 + a2 * S2 - a3 * S3
      dE3 <-  lambda3 * S3 - sigma * E3 + a2 * E2 - a3 * E3
      dI3 <-  sigma * E3 - gamma * I3 + a2 * I2 - a3 * I3
      dR3 <-  gamma * I3 - omega * R3 + a2 * R2 - a3 * R3

      dS4 <- - lambda4 * S4 + omega * R4 + a3 * S3 - a4 * S4
      dE4 <-  lambda4 * S4 - sigma * E4 + a3 * E3 - a4 * E4
      dI4 <-  sigma * E4 - gamma * I4 + a3 * I3 - a4 * I4
      dR4 <-  gamma * I4 - omega * R4 + a3 * R3 - a4 * I4  # NOTE: original had - a4 * I4 (possible bug: should be - a4 * R4)

      dS5 <- - lambda5 * S5 + omega * R5 + a4 * S4 - a5 * S5
      dE5 <-  lambda5 * S5 - sigma * E5 + a4 * E4 - a5 * E5
      dI5 <-  sigma * E5 - gamma * I5 + a4 * I4 - a5 * I5
      dR5 <-  gamma * I5 - omega * R5 + a4 * R4 - a5 * R5

      dS6 <- - lambda6 * S6 + omega * R6 + a5 * S5 - a6 * S6
      dE6 <-  lambda6 * S6 - sigma * E6 + a5 * E5 - a6 * E6
      dI6 <-  sigma * E6 - gamma * I6 + a5 * I5 - a6 * I6
      dR6 <-  gamma * I6 - omega * R6 + + a5 * R5 - a6 * R6

      dS7 <- - lambda7 * S7 + omega * R7 + a6 * S6 - a7 * S7
      dE7 <-  lambda7 * S7 - sigma * E7 + a6 * E6 - a7 * E7
      dI7 <-  sigma * E7 - gamma * I7 + a6 * I6 - a7 * I7
      dR7 <-  gamma * I7 - omega * R7 + a6 * R6 - a7 * R7

      dS8 <- - lambda8 * S8 + omega * R8 + a7 * S7 - a8 * S8
      dE8 <-  lambda8 * S8 - sigma * E8 + a7 * E7 - a8 * E8
      dI8 <-  sigma * E8 - gamma * I8 + a7 * I7 - a8 * I8
      dR8 <-  gamma * I8 - omega * R8 + a7 * R7 - a8 * R8

      dS9 <- - lambda9 * S9 + omega * R9 + a8 * S8 - mu_d * S9
      dE9 <-  lambda9 * S9 - sigma * E9 + a8 * E8 - mu_d * E9
      dI9 <-  sigma * E9 - gamma * I9 + a8 * I8 - mu_d * I9
      dR9 <-  gamma * I9 - omega * R9 + a8 * R8 - mu_d * R9 # no aging out and only deaths

      # incidence
      inc1 <- sigma * E1
      inc2 <- sigma * E2
      inc3 <- sigma * E3
      inc4 <- sigma * E4
      inc5 <- sigma * E5
      inc6 <- sigma * E6
      inc7 <- sigma * E7
      inc8 <- sigma * E8
      inc9 <- sigma * E9

      list(c(dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2, dS3, dE3, dI3, dR3,
             dS4, dE4, dI4, dR4, dS5, dE5, dI5, dR5, dS6, dE6, dI6, dR6,
             dS7, dE7, dI7, dR7, dS8, dE8, dI8, dR8, dS9, dE9, dI9, dR9),
           inc1 = inc1, inc2 = inc2, inc3 = inc3, inc4 = inc4, inc5 = inc5, inc6 = inc6, inc7 = inc7, inc8 = inc8, inc9 = inc9)
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
