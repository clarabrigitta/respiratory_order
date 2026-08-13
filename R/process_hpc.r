# combine the per-task HPC outputs (out1.rds ... out8.rds) into a single
# named list, matching what fit_model_rcpp.r used to save.

source(here::here("R", "library_hpc.r"))
source(here("R", "create_combinations.r"))

combinations <- create_combinations()

# folder the HPC job wrote to, i.e. the date the fit was run
date <- "12082026"

files <- here("inst", "outdata", date, paste0("out", seq_along(combinations), ".rds"))

results <- lapply(files, readRDS)
names(results) <- sapply(combinations, function(x) x$name)

saveRDS(results, file = here("inst", "outdata", paste0("parameters_", date)))
