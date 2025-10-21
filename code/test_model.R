# ---- Simple SEIR Model in R: facet grid over s (rows) and g (cols) ----

# Fixed parameters
N   <- 1000000
# Time (days)
times <- seq(0, 359, by = 1)

# Helper function: Define contact function over time (contacts per person per day)
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


# Helper function: converting day to month of year
day_to_month <- function(t) {
  d <- (floor(t) %% 365) + 1
  cum <- cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31))
  which(d <= cum)[1]
}

# SEIR ODEs
seir <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    cval <- c_t(t)           # contact rate at time t
    m <- month_foi[day_to_month(t)]  # month multiplier
    beta <- b * cval * m        # effective transmission rate
    lambda <- beta * I / N
    
    dS <- -lambda * S
    dE <-  lambda * S - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I
    
    Rt <- r0 * (S / N) * m
    
    list(c(dS, dE, dI, dR), 
         beta = beta, contacts = cval, month_mult = m, lambda = lambda, Rt = Rt)
  })
}

theta_to_mvec <- function(theta) {
  # identifiability: mean(m) = 1
  raw <- exp(theta - mean(theta))
  raw / mean(raw)
}
inv_logit <- function(x) 1/(1 + exp(-x))

# Grid of s and g
s_vals  <- 1:5
g_vals  <- 1:5
I0_vals <- c(1, 10, 20, 50, 100)
r0_vals <- c(0.5, 0.8, 1, 2, 3)

grid <- expand.grid(s = s_vals, g = g_vals,
                    I0 = I0_vals, r0 = r0_vals,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

# Run model for each (s, g)
sim_list <- pmap(grid, function(s, g, I0, r0) {
  sigma <- 1 / s          # rate from E to I
  gamma <- 1 / g         # rate from I to R
  b     <- r0 * gamma          
  Tg    <- 1/sigma + 1/gamma    # generation time
  
  # Initial conditions
  E0 <- 0
  R0 <- 0
  S0 <- N - I0 - E0 - R0
  y0 <- c(S = S0, E = E0, I = I0, R = R0)
  
  out <- ode(y = y0, times = times, func = seir,
             parms = list(N = N, sigma = sigma, gamma = gamma, b = b, r0 = r0))
  out <- as.data.frame(out)
  
  out %>% 
    mutate(s = s, g = g, Tg = Tg, I0 = I0, r0 = r0) %>% 
    pivot_longer(cols = c(S, E, I, R), names_to = "compartment", values_to = "count")
})

all_out <- bind_rows(sim_list)

# Plot a 5x5 grid: rows = s, cols = g
sub_out <- all_out %>% filter(s == 2, g == 5)
# sub_out <- all_out %>% filter(I0 == 10, r0 == 2)

ggplot(sub_out, aes(x = time, y = count, colour = compartment)) +
  geom_line(size = 0.4) +
  facet_grid(factor(I0, levels = I0_vals, labels = paste0("I0=", I0_vals)) ~
               factor(r0, levels = r0_vals, labels = paste0("r0=", r0_vals))) +
  # facet_grid(factor(s, levels = s_vals, labels = paste0("s=", s_vals)) ~
  #              factor(g, levels = g_vals, labels = paste0("g=", g_vals))) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size = 9)) +
  labs(
    x = "time (days)",
    y = "count",
    colour = "Compartment",
    subtitle = paste0("s = 2",
                      ", g = 5")
  )

# test incidence threshold
data <- sim_list[[422]] %>% filter(compartment == "I") %>% 
  mutate(epi_onset = rollapply(count, width = 5, FUN = function(v) all(diff(v) > 0), align = "right", fill = NA))

ggplot(data) +
  geom_line(aes(x = time, y = count)) +
  geom_point(data = subset(data, epi_onset == TRUE), aes(x = time, y = count), size = 2) +
  theme_bw(base_size = 10)

# Rt calculation
data <- sim_list[[422]] %>% 
  pivot_wider(names_from = compartment, values_from = count) %>%
  group_by(time) %>% 
  mutate(rt = r0 * (S/N))

ggplot(data) +
  geom_line(aes(x = time, y = rt)) +
  theme_bw(base_size = 10)
