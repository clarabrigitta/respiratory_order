# Parameter summary table for the fitted models.
# For each virus and parameter: MAP, 95% credible interval (equal-tailed),
# rank-normalised R-hat, bulk ESS and tail ESS (posterior package).
# Returns a long table with columns pathogen, parameter, metric, value.
#
# results : list of BayesianTools fits, one element per virus
# start   : first iteration kept (burn-in), as in process_outputs.r

process_parameters <- function(results, start = 3) {

  bind_rows(lapply(names(results), function(v) {
    # keep chains separate so R-hat/ESS can compare them
    chain <- getSample(results[[v]], start = start, parametersOnly = TRUE, coda = TRUE)
    map   <- MAP(results[[v]], start = start)$parametersMAP

    summarise_draws(as_draws_array(chain),
                    lower_95 = ~unname(quantile(.x, 0.025)),
                    upper_95 = ~unname(quantile(.x, 0.975)),
                    rhat     = posterior::rhat,
                    ess_bulk = posterior::ess_bulk,
                    ess_tail = posterior::ess_tail) %>%
      mutate(pathogen = v,
             MAP = unname(map[variable]))})) %>%
    rename(parameter = variable) %>%
    pivot_longer(c(MAP, lower_95, upper_95, rhat, ess_bulk, ess_tail),
                 names_to = "metric", values_to = "value") %>%
    select(pathogen, parameter, metric, value)
}

results <- readRDS(file = here("inst", "outdata", "parameters_18092026")) # change date as needed
params_summary <- process_parameters(results)
saveRDS(params_summary, file = here("inst", "outdata", "params_summary_04092026")) # change date as needed
