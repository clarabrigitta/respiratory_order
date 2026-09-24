# need to run/source: prepare_inputs.r, model_cpp.r and run fit_model_rcpp.r (line 1-107)

# Simulation-recovery: draw "true" parameter sets from the prior, simulate data
# from each, refit, and check how well the posteriors recover the truth.

age_groups <- c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")
param_names <- c(paste0("detection_", c("0to4","5to14","15to44","45to64","65plus")),
                 paste0("sus_",       c("0to4","5to14","15to44","45to64","65plus")),
                 "imm_duration", "log10_import", "p_inf")

n = 3 # RSV

# Prior ----
# used both to draw the true values and to refit
lb <- combinations[[n]]$lb
ub <- combinations[[n]]$ub
imm_meanlog <- log(combinations[[n]]$imm_period)
imm_sdlog <- combinations[[n]]$imm_sd

prior_local <- createPrior(
  density = function(param) {
    ld_uniform <- sum(dunif(param[u_idx], min = lb[u_idx], max = ub[u_idx], log = TRUE))
    ld_imm <- dlnorm(param[idx$imm],
                     meanlog = imm_meanlog, sdlog = imm_sdlog,
                     log = TRUE)
    ld_uniform + ld_imm
  },
  sampler = function(n = 1) {
    u <- mapply(function(lo, hi) runif(n, lo, hi), lb[u_idx], ub[u_idx])
    d <- rlnorm(n, meanlog = imm_meanlog, sdlog = imm_sdlog)
    # u columns follow u_idx: det(1-5), sus(6-10), log10_import, p_inf
    if (n == 1) c(u[1:10], d, u[11:12]) else cbind(u[, 1:10], d, u[, 11:12])
  },
  lower = c(lb[1:10], 0,   lb[12:13]),
  upper = c(ub[1:10], Inf, ub[12:13])
)

# Draw true parameter sets ----
n_sims <- 100

set.seed(24)
theta_draws <- prior_local$sampler(n_sims)
colnames(theta_draws) <- param_names

# Simulate data ----
simulate_data <- function(theta) {
  # true infections
  infections <- run_model_rcpp(p_sus_bands = theta[idx$sus],
                               n = n,
                               imm_days = theta[idx$imm],
                               imports = 10^theta[idx$imp],
                               p_inf   = theta[idx$pinf])

  # add detection + noise
  det_infections <- infections * rep(theta[idx$det], each = nrow(infections))
  sim_obs <- matrix(rpois(length(det_infections), det_infections), nrow = nrow(infections))
  dimnames(sim_obs) <- list(rownames(infections), age_groups)
  sim_obs
}

plot_sim <- function(sim_obs) {
  as.data.frame(sim_obs) %>%
    mutate(date = as.Date(rownames(sim_obs))) %>%
    pivot_longer(cols = -date, names_to = "age_group", values_to = "Count") %>%
    mutate(Pathogen = "RSV") %>%
    select(Pathogen, date, age_group, Count) %>%
    ggplot()+
    geom_line(aes(x = date, y = Count, colour = age_group)) +
    scale_color_viridis_d(option = "H") +
    scale_x_date(date_breaks = "3 month") +
    labs(x = "Weeks", colour = "Age group") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11))
}

plot_all_sims <- function(sim_obs_list) {
  bind_rows(lapply(seq_along(sim_obs_list), function(k) {
    as.data.frame(sim_obs_list[[k]]) %>%
      mutate(sim = k, date = as.Date(rownames(sim_obs_list[[k]])))
  })) %>%
    pivot_longer(cols = all_of(age_groups), names_to = "age_group", values_to = "Count") %>%
    mutate(age_group = factor(age_group, levels = age_groups)) %>%
    ggplot()+
    geom_line(aes(x = date, y = Count, group = sim), alpha = 0.3) +
    facet_wrap(~age_group, scales = "free_y") +
    scale_x_date(date_breaks = "3 month") +
    labs(x = "Weeks") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# simulate every dataset up front (cheap), so they can be viewed before fitting
sim_obs_list <- lapply(seq_len(n_sims), function(k) {
  set.seed(24 + k)
  simulate_data(theta_draws[k, ])
})

plot_all_sims(sim_obs_list)  # all datasets
plot_sim(sim_obs_list[[1]])  # single dataset
theta_draws[1, ]             # its true parameters

# Fit to simulated dataset ---- (only needed for recovery; skip to just view data)
settings <- list(iterations = 500000, burnin = 200000, nrChains = 1)

fit_sim <- function(sim_obs) {
  obs_vec <- as.vector(sim_obs)

  likelihood <- function(param) {
    detection_rates <- param[idx$det]
    p_sus_bands     <- param[idx$sus]
    imm_duration    <- param[idx$imm]
    imports         <- 10^param[idx$imp]
    p_inf           <- param[idx$pinf]

    model_out <- run_model_rcpp(p_sus_bands, n, imm_duration, imports, p_inf)

    expected_vec <- as.vector(model_out) *
      rep(detection_rates, each = nrow(model_out))

    ll <- sum(dpois(obs_vec, expected_vec, log = TRUE))
    if (!is.finite(ll)) return(-1e10)
    ll
  }

  setup <- createBayesianSetup(
    likelihood = likelihood,
    prior = prior_local,
    parallel = FALSE,
    names = param_names
  )

  chains <- mclapply(1:4,
                     function(x) {
                       runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)
                     },
                     mc.cores = 4)
  createMcmcSamplerList(chains)
}

# fit one replicate; k indexes theta_draws and sim_obs_list
run_one_sim <- function(k) {
  set.seed(24 + k)
  theta_true <- theta_draws[k, ]
  sim_obs <- sim_obs_list[[k]]
  out <- fit_sim(sim_obs)

  saveRDS(list(theta_true = theta_true, sim_obs = sim_obs, out = out),
          file = here("inst", "outdata", paste0("sim_", combinations[[n]]$name, "_", k)))

  # keep only a thinned posterior in memory (full chains are on disk).
  # start = 3 drops each chain's start value: it is a prior draw, and post_sd
  # below is not robust to it (sd inflated up to ~200x for tight parameters).
  list(k = k, theta_true = theta_true, sim_obs = sim_obs,
       post = getSample(out, start = 3, thin = 10))
}

# on HPC, run one replicate per array task instead, e.g.
# run_one_sim(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
sims <- lapply(seq_len(n_sims), function(k) {
  cat("Simulation", k, "of", n_sims, "-", format(Sys.time(), "%H:%M:%S"), "\n")
  run_one_sim(k)
})

# Assess recovery ----
prior_sd <- setNames(apply(prior_local$sampler(10000), 2, sd), param_names)

recovery <- bind_rows(lapply(sims, function(s) {
  tibble(sim       = s$k,
         parameter = factor(param_names, levels = param_names),
         true      = s$theta_true,
         median    = apply(s$post, 2, median),
         lower     = apply(s$post, 2, quantile, probs = 0.025),
         upper     = apply(s$post, 2, quantile, probs = 0.975),
         post_sd   = apply(s$post, 2, sd),
         # fraction of posterior draws below the truth: ~Uniform(0, 1) across sims if calibrated
         rank      = colMeans(sweep(s$post, 2, s$theta_true, "<")))
}))

recovery %>%
  mutate(prior_sd = prior_sd[as.character(parameter)]) %>%
  group_by(parameter) %>%
  summarise(coverage_95 = mean(true >= lower & true <= upper), # should be ~0.95
            mean_z      = mean((median - true) / post_sd),     # bias, should be ~0
            contraction = mean(1 - post_sd^2 / prior_sd^2))     # ~1 = data informative, ~0 = prior only

# true vs estimated
ggplot(recovery, aes(x = true, y = median)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2) +
  facet_wrap(~parameter, scales = "free") +
  labs(x = "True value", y = "Posterior median (95% CrI)") +
  theme_classic()

# calibration: histograms should look roughly flat
ggplot(recovery, aes(x = rank)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1)) +
  facet_wrap(~parameter) +
  labs(x = "Fraction of posterior draws below true value", y = "Simulations") +
  theme_classic()
