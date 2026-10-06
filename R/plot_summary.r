# Summary figures for the fitted models.
# Each figure lives in its own script under R/plots/; this script loads the
# shared inputs and calls them.
#
source(here::here("R", "library_hpc.r"))
source(here("R", "create_combinations.r"))

combinations <- create_combinations()

# scenario to process: array task i -> i-th scenario (first one if run locally)
scenario_names <- sapply(combinations, function(x) x$scenario)
scenario <- unique(scenario_names)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]
combinations <- combinations[scenario_names == scenario]

# define date for file use/saving: the run date set in inst/bash/run.sh (today if run without it)
date <- Sys.getenv("RUN_DATE", format(Sys.Date(), "%d%m%Y"))
plot_dir <- here("inst", "plots", date, scenario)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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

# trajectory and Rt plots ----
# these only need the small traj/rt files, so they run before the large MCMC
# output is loaded
results_traj <- readRDS(file = here("inst", "outdata", date, paste0("traj_", date, "_", scenario)))

# panel plot of trajectories for all pathogens, needs the NON-age-stratified `combined`
plot_traj_total(results_traj, data, combined = mean_total_scotland, age_groups, pathogen_map, outdir = plot_dir, date = date)

# age-specific, needs the AGE-STRATIFIED `combined`
plot_traj_age(results_traj, data, combined = mean_age_scotland, age_groups, pathogen_map, outdir = plot_dir, date = date)

rm(results_traj)

results_rt <- readRDS(file = here("inst", "outdata", date, paste0("rt_", date, "_", scenario)))
plot_rt(results_rt, outdir = plot_dir, date = date)
rm(results_rt)

# parameter plots, one saved figure of each type per virus ----
results <- readRDS(file = here("inst", "outdata", date, paste0("parameters_", date, "_", scenario)))

for (virus_name in names(results)) {
  ggsave(filename = file.path(plot_dir, paste0("prior_posterior_", virus_name, "_", date, ".png")),
         plot = plot_prior_posterior(results, combinations, virus_name), width = 15, height = 9, dpi = 300)

  ggsave(filename = file.path(plot_dir, paste0("mcmc_trace_", virus_name, "_", date, ".png")),
         plot = plot_mcmc_trace(results, virus_name), width = 15, height = 9, dpi = 300)

  ggsave(filename = file.path(plot_dir, paste0("correlation_", virus_name, "_", date, ".png")),
         plot = plot_correlation(results, virus_name), width = 16, height = 16, dpi = 300)

  gc()
}
