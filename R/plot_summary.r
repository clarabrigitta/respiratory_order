# Summary figures for the fitted models.
# Each figure lives in its own script under R/plots/; this script loads the
# shared inputs and calls them.
#
source(here::here("R", "library_hpc.r"))
source(here("R", "create_combinations.r"))

combinations <- create_combinations()

# define date for file use/saving
date <- format(Sys.Date(), "%d%m%Y")
dir.create(here("inst", "plots", date), recursive = TRUE, showWarnings = FALSE)

# source relevant data ----
pathogen_map <- c(fluA = "Influenza A",
                  fluB = "Influenza B",
                  RSV  = "RSV",
                  hCOV = "Seasonal coronavirus",
                  AdV  = "Adenovirus",
                  RV   = "Rhinovirus",
                  hMPV = "HMPV",
                  PIV  = "Parainfluenza (Any Type)")

# load mean fortnightly number of contacts in Scotland (total and age-breakdown)
inputs <- readRDS(here("inst", "outdata", "hpc_inputs.rds"))
mean_total_scotland <- inputs$mean_total_scotland
mean_age_scotland <- inputs$mean_age_scotland

# load plotting functions ----
for (f in list.files(here("R", "plots"), pattern = "[.][rR]$", full.names = TRUE)) source(f)

# load model output ----
results <- readRDS(file = here("inst", "outdata", paste0("parameters_", date)))
results_traj <- readRDS(file = here("inst", "outdata", paste0("traj_", date)))
results_rt <- readRDS(file = here("inst", "outdata", paste0("rt_", date)))

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

# fitted parameter names
params <- c(paste0("detection_", c("0to4","5to14","15to44","45to64","65plus")),
            paste0("sus_",       c("0to4","5to14","15to44","45to64","65plus")),
            "imm_duration", "log10_import", "p_inf")

# panel plot of trajectories for all pathogens ----
# needs the NON-age-stratified `combined`
fig_traj <- plot_traj_total(results_traj, data, combined = mean_total_scotland, age_groups, pathogen_map)

# age-specific ----
# needs the AGE-STRATIFIED `combined`
figs_traj_age <- plot_traj_age(results_traj, data, combined = mean_age_scotland, age_groups, pathogen_map)

# prior vs posterior of each parameter, one saved figure per virus ----
# needs `combinations` (priors are rebuilt from it; the fit's own prior closure
# does not survive saveRDS)
figs_prior_posterior <- list()

for (virus_name in names(results)) {
  fig_prior_posterior <- plot_prior_posterior(results, combinations, virus_name)

  ggsave(filename = here("inst", "plots", date, paste0("prior_posterior_", virus_name, "_", date, ".png")),
         plot = fig_prior_posterior, width = 15, height = 9, dpi = 300)

  figs_prior_posterior[[virus_name]] <- fig_prior_posterior
}

# MCMC traceplots of each parameter, one saved figure per virus ----
figs_mcmc_trace <- list()

for (virus_name in names(results)) {
  fig_mcmc_trace <- plot_mcmc_trace(results, virus_name)

  ggsave(filename = here("inst", "plots", date, paste0("mcmc_trace_", virus_name, "_", date, ".png")),
         plot = fig_mcmc_trace, width = 15, height = 9, dpi = 300)

  figs_mcmc_trace[[virus_name]] <- fig_mcmc_trace
}

# pairwise parameter correlations, one saved figure per virus ----
figs_correlation <- list()

for (virus_name in names(results)) {
  fig_correlation <- plot_correlation(results, virus_name)

  ggsave(filename = here("inst", "plots", date, paste0("correlation_", virus_name, "_", date, ".png")),
         plot = fig_correlation, width = 16, height = 16, dpi = 300)

  figs_correlation[[virus_name]] <- fig_correlation
}

# plot Rt trajectories ----
fig_rt <- plot_rt(results_rt)
