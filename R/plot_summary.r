# Summary figures for the fitted models.
# Each figure lives in its own script under R/plots/; this script loads the
# shared inputs and calls them.
#
# Assumes already sourced/available: library_and_scripts.r, model_rcpp.r,
# fit_model_rcpp.r (pathogen_map, times, tots, daily_births, hazard_death,
# contacts_prepped, combinations) and `combined` from explore_contacts.r.

# load plotting functions ----
for (f in list.files(here("R", "plots"), pattern = "[.][rR]$", full.names = TRUE)) source(f)

# load model output ----
results <- readRDS(file = here("inst", "outdata", "parameters_18092026"))
results_traj <- readRDS(file = here("inst", "outdata", "traj_18092026"))

# load case data ----
data <- read_csv("inst/data/cases_all_respiratory_pathogens_by_agegroup_sex_20251008.csv") %>%
  filter(Sex == "Total", !AgeGroup %in% c("Total", "Unknown"), !Pathogen %in% c("Influenza (All)", "COVID-19", "Mycoplasma pneumoniae")) %>%
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>%
  filter(WeekBeginning >= as.Date("2020-03-23"),
         WeekBeginning <= as.Date("2022-03-02")) %>%
  pivot_wider(names_from = AgeGroup, values_from = NumberCasesPerWeek) %>%
  mutate(`0 to 4` = `<1` + `1 to 4`,
         `65+` = `65 to 74` + `75+`) %>%
  select(-c(`<1`, `1 to 4`, `65 to 74`, `75+`)) %>%
  arrange(WeekBeginning)

# load age group bands
age_groups  <- c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")
# age_colours <- setNames(viridis(5), age_groups)
# need to also load combined from explore_contacts.r either total or age-stratified (lines 271 or 295)

# fitted parameter names
params <- c(paste0("detection_", c("0to4","5to14","15to44","45to64","65plus")),
            paste0("sus_",       c("0to4","5to14","15to44","45to64","65plus")),
            "imm_duration", "log10_import", "p_inf")

# panel plot of trajectories for all pathogens ----
# needs the NON-age-stratified `combined` (explore_contacts.r line 271)
fig_traj <- plot_traj_total(results_traj, data, combined, age_groups, pathogen_map)

# age-specific ----
# needs the AGE-STRATIFIED `combined` (explore_contacts.r line 295)
figs_traj_age <- plot_traj_age(results_traj, data, combined, age_groups, pathogen_map)

# MCMC traceplot panel for all viruses and parameters ----
plot_chains(results, params)

# panel plot of posterior distribution of parameter estimates ----
fig_postdist <- plot_posterior_distribution(results)

# prior vs posterior of each parameter, one saved figure per virus ----
# needs `combinations` (priors are rebuilt from it; the fit's own prior closure
# does not survive saveRDS)
figs_prior_posterior <- list()

for (virus_name in names(results)) {
  fig_prior_posterior <- plot_prior_posterior(results, combinations, virus_name)

  ggsave(filename = here("inst", "plots", paste0("prior_posterior_", virus_name, "_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_prior_posterior, width = 15, height = 9, dpi = 300)

  figs_prior_posterior[[virus_name]] <- fig_prior_posterior
}

# MCMC traceplots of each parameter, one saved figure per virus ----
figs_mcmc_trace <- list()

for (virus_name in names(results)) {
  fig_mcmc_trace <- plot_mcmc_trace(results, virus_name)

  ggsave(filename = here("inst", "plots", paste0("mcmc_trace_", virus_name, "_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_mcmc_trace, width = 15, height = 9, dpi = 300)

  figs_mcmc_trace[[virus_name]] <- fig_mcmc_trace
}

# pairwise parameter correlations, one saved figure per virus ----
figs_correlation <- list()

for (virus_name in names(results)) {
  fig_correlation <- plot_correlation(results, virus_name)

  ggsave(filename = here("inst", "plots", paste0("correlation_", virus_name, "_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_correlation, width = 16, height = 16, dpi = 300)

  figs_correlation[[virus_name]] <- fig_correlation
}

# panel plot of range of parameter estimates ----
fig_paramsbox <- plot_parameters_boxplot(results, params)

# plot Rt trajectories ----
fig_rt <- plot_rt(results_rt)

get_draws <- function(results, start = 3, thin = "auto") {
  bind_rows(lapply(names(results), function(v) {
    getSample(results[[v]], start = start, thin = thin, parametersOnly = TRUE) %>%
      as.data.frame() %>% mutate(pathogen = v)
  })) %>%
    pivot_longer(-pathogen, names_to = "parameter", values_to = "value") %>%
    # imports are stored as log10 already; everything else is positive, so log it
    mutate(parameter = if_else(parameter == "log10_import", "import", parameter),
           value = if_else(parameter == "import", value,
                           log10(if_else(value > 0, value, NA_real_)))) %>%
    filter(is.finite(value))
}

coda_density <- function(y, n = 512) {
  y  <- y[!is.na(y)]
  bw <- 1.06 * min(sd(y), IQR(y) / 1.34) * length(y)^-0.2
  d  <- density(y, width = 4 * bw, n = n)   # log scale: no boundary to reflect at
  data.frame(x = d$x, y = d$y)
}

plot_posterior_density <- function(results, start = 3, thin = "auto", normalise = TRUE) {
  dens <- get_draws(results, start, thin) %>%
    group_by(pathogen, parameter) %>%
    reframe(coda_density(value)) %>%           # dplyr >= 1.1
    group_by(pathogen, parameter) %>%
    mutate(y = if (normalise) y / max(y) else y) %>%
    ungroup()
  
  ggplot(dens, aes(x = x, y = y, colour = pathogen, fill = pathogen)) +
    geom_area(alpha = 0.15, position = "identity", linewidth = 0) +
    geom_line(linewidth = 0.7) +
    scale_colour_viridis_d() + scale_fill_viridis_d() +
    scale_x_continuous(n.breaks = 5) +
    facet_wrap(~parameter, scales = "free", nrow = 2) +
    labs(x = expression(log[10]~"(value)"),
         y = if (normalise) "density (scaled to peak)" else "density") +
    theme_classic() +
    theme(legend.title = element_blank(), strip.text = element_text(face = "bold"))
}
plot_posterior_density(results)


# plot susceptible S(t) by age group for each virus, with mean contacts on 2nd axis ----
# needs the AGE-STRATIFIED `combined` (explore_contacts.r line 295)
fig_susceptible <- plot_susceptible_age(results, combined, combinations, age_groups,
                                        times, tots, daily_births, hazard_death, contacts_prepped)

# plot TOTAL susceptible S(t) for each virus, with mean contacts on 2nd axis ----
# needs the NON-age-stratified `combined` (explore_contacts.r line 271)
fig_susceptible_total <- plot_susceptible_total(results, combined, combinations,
                                                times, tots, daily_births, hazard_death, contacts_prepped)