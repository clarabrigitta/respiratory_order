# age-stratified SEIRS model ----

## helper data and functions ----
### date data frame for assistance/reference
dates <- data.frame(date = seq(as.Date("23-03-2020", format = "%d-%m-%Y"), as.Date("02-03-2022", format = "%d-%m-%Y"), 1)) %>% 
  mutate(time = 0:(n()-1),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date)) %>% 
  mutate(fortnight_n = as.integer(factor(fortnight)))

### define function to choose which fortnightly contact matrix to use (choose the contact matrix based on a selection from fortnight_matrix)
c_t <- function(t){
  #filter fortnight_n value/row based on t, t defined in ode but is continuous and not discrete, using findInterval to map to discrete
  fortnight_value <- filter(dates, time == findInterval(t, dates$time)) %>%
    select(fortnight_n) %>%
    as.integer()
  
  #return correct fortnight_matrix contact matrix based on fortnight_n value
  return(fortnight_matrix[[fortnight_value]]$matrix)
}

### create data frame defining model parameters
source(here("R", "create_combinations.R"))
combinations <- create_combinations()

### generate functional form of immunity given age (distribution) to use for how many people in R compartment for each age group
imm_distribution <- function(a, mu, sigma){
  exp(- (a - mu)^2 / (2 * sigma^2))
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
      dS1 <- - lambda1 * S1 + omega * R1
      dE1 <-  lambda1 * S1 - sigma * E1
      dI1 <-  sigma * E1 - gamma * I1
      dR1 <-  gamma * I1 - omega * R1

      dS2 <- - lambda2 * S2 + omega * R2
      dE2 <-  lambda2 * S2 - sigma * E2
      dI2 <-  sigma * E2 - gamma * I2
      dR2 <-  gamma * I2 - omega * R2
     
      dS3 <- - lambda3 * S3 + omega * R3
      dE3 <-  lambda3 * S3 - sigma * E3
      dI3 <-  sigma * E3 - gamma * I3
      dR3 <-  gamma * I3 - omega * R3

      dS4 <- - lambda4 * S4 + omega * R4
      dE4 <-  lambda4 * S4 - sigma * E4
      dI4 <-  sigma * E4 - gamma * I4
      dR4 <-  gamma * I4 - omega * R4

      dS5 <- - lambda5 * S5 + omega * R5
      dE5 <-  lambda5 * S5 - sigma * E5
      dI5 <-  sigma * E5 - gamma * I5
      dR5 <-  gamma * I5 - omega * R5

      dS6 <- - lambda6 * S6 + omega * R6
      dE6 <-  lambda6 * S6 - sigma * E6
      dI6 <-  sigma * E6 - gamma * I6
      dR6 <-  gamma * I6 - omega * R6

      dS7 <- - lambda7 * S7 + omega * R7
      dE7 <-  lambda7 * S7 - sigma * E7
      dI7 <-  sigma * E7 - gamma * I7
      dR7 <-  gamma * I7 - omega * R7

      dS8 <- - lambda8 * S8 + omega * R8
      dE8 <-  lambda8 * S8 - sigma * E8
      dI8 <-  sigma * E8 - gamma * I8
      dR8 <-  gamma * I8 - omega * R8

      dS9 <- - lambda9 * S9 + omega * R9
      dE9 <-  lambda9 * S9 - sigma * E9
      dI9 <-  sigma * E9 - gamma * I9
      dR9 <-  gamma * I9 - omega * R9
      
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

## run model ----

### average population estimate between 2020-2022 for each age group
tot1 <- 252675
tot2 <- 405926
tot3 <- 350290
tot4 <- 801593
tot5 <- 699525
tot6 <- 667172
tot7 <- 800879
tot8 <- 676289
tot9 <- 768583

### starting compartment split for each age group
y0 <- c(S1 = tot1-5-floor(tot1*(1- imm_distribution(1, 3, 10))), E1 = 0, I1 = 5, R1 = floor(tot1*(1- imm_distribution(1, 3, 10))), # using random normal distribution for now, not informed by anything (need to change into something proper later)
        S2 = tot2-5-floor(tot2*(1- imm_distribution(2, 3, 10))), E2 = 0, I2 = 5, R2 = floor(tot2*(1- imm_distribution(2, 3, 10))),
        S3 = tot3-5-floor(tot3*(1- imm_distribution(3, 3, 10))), E3 = 0, I3 = 5, R3 = floor(tot3*(1- imm_distribution(3, 3, 10))),
        S4 = tot4-5-floor(tot4*(1- imm_distribution(4, 3, 10))), E4 = 0, I4 = 5, R4 = floor(tot4*(1- imm_distribution(4, 3, 10))),
        S5 = tot5-5-floor(tot5*(1- imm_distribution(5, 3, 10))), E5 = 0, I5 = 5, R5 = floor(tot5*(1- imm_distribution(5, 3, 10))),
        S6 = tot6-5-floor(tot6*(1- imm_distribution(6, 3, 10))), E6 = 0, I6 = 5, R6 = floor(tot6*(1- imm_distribution(6, 3, 10))),
        S7 = tot7-5-floor(tot7*(1- imm_distribution(7, 3, 10))), E7 = 0, I7 = 5, R7 = floor(tot7*(1- imm_distribution(7, 3, 10))),
        S8 = tot8-5-floor(tot8*(1- imm_distribution(8, 3, 10))), E8 = 0, I8 = 5, R8 = floor(tot8*(1- imm_distribution(8, 3, 10))),
        S9 = tot9-5-floor(tot9*(1- imm_distribution(9, 3, 10))), E9 = 0, I9 = 5, R9 = floor(tot9*(1- imm_distribution(9, 3, 10))))

### duration of model run
times <- seq(0, 705, by = 1)

### parallel run model
print(paste("start:", Sys.time()))
out <- mclapply(1:7,
                function(n) {
                  ode(y = y0, times = times, func = seirs, 
                      parms = list(sigma = 1 / combinations[[n]]$inc_period, gamma = 1 / combinations[[n]]$inf_period, omega = 1 / combinations[[n]]$imm_period, p_inf = combinations[[n]]$p_inf))
                },
                mc.cores = 4)
print(paste("end:", Sys.time()))

### name list output
names(out) <- map_chr(combinations, 1)
                
## plot model output ----
fig <- out %>% 
  as.data.frame() %>% 
  rename(inc1 = inc1.E1, inc2 = inc2.E2, inc3 = inc3.E3, inc4 = inc4.E4, inc5 = inc5.E5, inc6 = inc6.E6, inc7 = inc7.E7, inc8 = inc8.E8, inc9 = inc9.E9) %>%
  pivot_longer(cols = 2:ncol(out), names_to = "compartment", values_to = "count") %>% 
  # mutate(compartment = factor(compartment, levels = c("S1", "E1", "I1", "R1", "inc1", "S2", "E2", "I2", "R2", "inc2"))) %>%
  # filter(compartment %in% c("S1", "E1", "I1", "R1", "inc1")) %>%
  filter(grepl('inc', compartment)) %>% 
  left_join(dates, by = join_by(time)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = count, colour = compartment)) +
  # scale_colour_manual(values = c(paletteer_c("ggthemes::Blue", 5), paletteer_c("ggthemes::Red", 5))) +
  scale_colour_viridis_d(option = "G", 
                         end = 0.8,
                         labels = c("inc1" = "0-4",
                                   "inc2" = "5-11",
                                   "inc3" = "12-17",
                                   "inc4" = "18-29",
                                   "inc5" = "30-39",
                                   "inc6" = "40-49",
                                   "inc7" = "50-59",
                                   "inc8" = "60-69",
                                   "inc9" = "70+")) +
  scale_x_date(date_breaks = "1 month") +
  labs(colour = "Age Group", y = "Count", x = "Time") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

## save plot
ggsave(filename = here("inst", "plots", "contact_fortnight_hMPV_130126.png"), plot = fig, width = 13, height = 9, dpi = 300)

## aggregate incidence from output
out %>% 
  as.data.frame() %>% 
  select(matches("inc|time")) %>% 
  pivot_longer(cols = 2:10) %>% 
  group_by(time) %>% 
  summarise(total_inc = sum(value)) %>% 
  left_join(dates, by = join_by(time)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = total_inc)) +
  scale_x_date(date_breaks = "1 month") +
  labs(y = "Count", x = "Time") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))
