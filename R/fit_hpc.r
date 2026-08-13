source(here::here("R", "library_hpc.r"))

n <- Sys.getenv("SLURM_ARRAY_TASK_ID")
n <- as.integer(n)

# fixed inputs built by R/prepare_inputs.r
inputs <- readRDS(here("inst", "outdata", "hpc_inputs.rds"))
fortnight_matrix <- inputs$fortnight_matrix
scot_population  <- inputs$scot_population
scot_births      <- inputs$scot_births

source(here("R", "model_rcpp.r"))

# helpers ----

fit_start <- as.Date("2020-03-23")
fit_end   <- as.Date("2022-03-02")
times <- seq(0, as.integer(fit_end - fit_start), by = 1)

contacts_prepped <- prepare_contacts_cpp(fortnight_matrix, fortnight_lookup)

pathogen_map <- c(fluA = "Influenza A",
                  fluB = "Influenza B",
                  RSV  = "RSV",
                  hCOV = "Seasonal coronavirus",
                  AdV  = "Adenovirus",
                  RV   = "Rhinovirus",
                  hMPV = "HMPV",
                  PIV  = "Parainfluenza (Any Type)")

virus_name    <- combinations[[n]]$name
pathogen_name <- pathogen_map[[virus_name]]

data <- read_csv(here("inst", "data", "cases_all_respiratory_pathogens_by_agegroup_sex_20251008.csv")) %>%
  filter(Sex == "Total",
         Pathogen == pathogen_name,
         !AgeGroup %in% c("Total", "Unknown")) %>%
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>%
  filter(WeekBeginning >= fit_start,
         WeekBeginning <= fit_end) %>%
  pivot_wider(names_from = AgeGroup, values_from = NumberCasesPerWeek) %>%
  mutate(`0 to 4` = `<1` + `1 to 4`,
         `65+` = `65 to 74` + `75+`) %>%
  select(-c(`<1`, `1 to 4`, `65 to 74`, `75+`)) %>%
  arrange(WeekBeginning) %>%
  select(Pathogen, "0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")

obs_vec <- as.vector(as.matrix(data[, -1]))

dates_seq  <- seq.Date(fit_start, fit_end, by = 1)
week_group <- as.character(floor_date(dates_seq, unit = "week", week_start = 1))

tots <- c(tot1, tot2, tot3, tot4, tot5, tot6, tot7, tot8, tot9)

# indexing for parameters (to keep better track)
idx <- list(det = 1:5, sus = 6:10, imm = 11, imp = 12)
band_of_group <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)   # 9 model groups -> 5 data bands
u_idx <- c(idx$det, idx$sus, idx$imp)
inc_cols <- 1 + 4 * 9 + seq_len(9)

# model runner ----

run_model_rcpp <- function(p_sus_bands, n, imm_days, imports = 1) {
  sigma <- 1 / combinations[[n]]$inc_period
  gamma <- 1 / combinations[[n]]$inf_period
  bg    <- imports / sum(tots)

  R_init <- tots * (1 - p_sus_bands[band_of_group])
  S_init <- tots - R_init
  E_init <- bg * S_init / sigma
  I_init <- bg * S_init / gamma

  y0_mat <- rbind(S = S_init - E_init - I_init, E = E_init,
                  I = I_init, R = R_init)
  y0 <- setNames(as.vector(y0_mat),
                 paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4)))

  out <- seirs_rcpp(y0 = y0,
                    times = times,
                    parms = list(sigma = sigma,
                                 gamma = gamma,
                                 omega = 1 / imm_days,
                                 p_inf = combinations[[n]]$p_inf,
                                 bg    = bg,
                                 mu_b  = daily_births,
                                 mu_d  = hazard_death),
                    contacts_prepped = contacts_prepped)

  daily_inc <- out[, inc_cols, drop = FALSE]

  daily_agg <- cbind(daily_inc[, 1],
                     daily_inc[, 2] + daily_inc[, 3],
                     daily_inc[, 4] + daily_inc[, 5] + daily_inc[, 6],
                     daily_inc[, 7] + daily_inc[, 8],
                     daily_inc[, 9])

  rowsum(daily_agg, week_group)
}

# MCMC fitting ----

lb <- combinations[[n]]$lb
ub <- combinations[[n]]$ub
imm_meanlog <- log(combinations[[n]]$imm_period)
imm_sdlog   <- combinations[[n]]$imm_sd

prior <- createPrior(
  density = function(param) {
    ld_uniform <- sum(dunif(param[u_idx], min = lb[u_idx], max = ub[u_idx], log = TRUE))
    ld_imm <- dlnorm(param[idx$imm], meanlog = imm_meanlog, sdlog = imm_sdlog, log = TRUE)
    ld_uniform + ld_imm
  },
  sampler = function(n = 1) {
    u <- mapply(function(lo, hi) runif(n, lo, hi), lb[u_idx], ub[u_idx])
    d <- rlnorm(n, meanlog = imm_meanlog, sdlog = imm_sdlog)
    if (n == 1) c(u[1:10], d, u[11]) else cbind(u[, 1:10], d, u[, 11])
  },
  lower = c(lb[1:10], 0,   lb[12]),
  upper = c(ub[1:10], Inf, ub[12])
)

# poisson likelihood function
likelihood <- function(param) {
  detection_rates <- param[idx$det]
  p_sus_bands     <- param[idx$sus]
  imm_duration    <- param[idx$imm]
  imports         <- 10^param[idx$imp]

  model_out <- run_model_rcpp(p_sus_bands, n, imm_duration, imports)
  expected_vec <- as.vector(model_out) * rep(detection_rates, each = nrow(model_out))

  ll <- sum(dpois(obs_vec, expected_vec, log = TRUE))
  if (!is.finite(ll)) return(-1e10)
  ll
}

# fitting setup
setup <- createBayesianSetup(likelihood = likelihood,
                             prior = prior,
                             parallel = FALSE,
                             names = c(paste0("detection_", c("0to4", "5to14", "15to44", "45to64", "65plus")),
                                       paste0("sus_", c("0to4", "5to14", "15to44", "45to64", "65plus")),
                                       "imm_duration", "log10_import"))

settings <- list(iterations = 500000, nrChains = 1, message = TRUE, burnin = 200000)

print(paste("start iteration number", n, virus_name, "time:", Sys.time()))

# fit model and save output
set.seed(24)
results <- mclapply(1:4,
                    function(x) {
                      runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)
                    },
                    mc.cores = 4)

out <- createMcmcSamplerList(results)

dir.create(here("inst", "outdata", format(Sys.Date(), "%d%m%Y")), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, file = here("inst", "outdata", format(Sys.Date(), "%d%m%Y"), paste0("out", n, ".rds")))

print(paste("end iteration number", n, virus_name, "time:", Sys.time()))
