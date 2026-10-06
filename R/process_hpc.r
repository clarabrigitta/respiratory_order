# combine the per-task HPC outputs (out1.rds ... out8.rds) into a single
# named list, matching what fit_model_rcpp.r used to save.

source(here::here("R", "library_hpc.r"))
source(here("R", "create_combinations.r"))

combinations <- create_combinations()

# scenario to process: array task i -> i-th scenario (first one if run locally)
scenario_names <- sapply(combinations, function(x) x$scenario)
scenario <- unique(scenario_names)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]

# folder the HPC job wrote to: the run date set in inst/bash/run.sh (today if run without it)
date <- Sys.getenv("RUN_DATE", format(Sys.Date(), "%d%m%Y"))

files <- here("inst", "outdata", date, paste0("out", which(scenario_names == scenario), ".rds"))

results <- lapply(files, readRDS)
names(results) <- sapply(combinations[scenario_names == scenario], function(x) x$name)

saveRDS(results, file = here("inst", "outdata", date, paste0("parameters_", date, "_", scenario)))
