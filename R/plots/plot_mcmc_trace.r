# MCMC traceplots for one virus, one facet per parameter, one line per chain.
# Mirrors the trace half of BayesianTools' plot() (coda::traceplot): raw values
# against iteration, chains kept separate. Returns the ggplot object.
#
#   plot_mcmc_trace(results, "RSV")
#
# results : list of BayesianTools fits, one element per virus
# start   : first iteration kept (burn-in), as in process_parameters.r
# thin    : "auto" keeps ~5000 draws per chain, as plot() does

plot_mcmc_trace <- function(results, virus, start = 3, thin = "auto", ncol = 5) {

  stopifnot(virus %in% names(results))

  chains <- getSample(results[[virus]], start = start, thin = thin,
                      parametersOnly = TRUE, coda = TRUE)
  par_names <- colnames(chains[[1]])

  data <- bind_rows(lapply(seq_along(chains), function(i) {
    # time() errors on these chains: BayesianTools stores an `end` that is not a
    # whole number of thinning steps after `start`, so build the index directly
    mcpar <- attr(chains[[i]], "mcpar")
    as.data.frame(as.matrix(chains[[i]])) %>%
      mutate(iteration = mcpar[1] + mcpar[3] * (seq_len(n()) - 1),
             chain = i)})) %>%
    pivot_longer(-c(iteration, chain), names_to = "parameter", values_to = "value") %>%
    mutate(parameter = factor(parameter, levels = par_names),
           chain = factor(chain))

  ggplot(data, aes(x = iteration, y = value, colour = chain)) +
    geom_line(linewidth = 0.2, alpha = 0.8) +
    scale_colour_viridis_d() +
    scale_x_continuous(n.breaks = 4, labels = scales::label_comma()) +
    labs(title = virus, x = "Iteration", y = "Value", colour = "Chain") +
    facet_wrap(~ parameter, scales = "free_y", ncol = ncol) +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "none",
          plot.title = element_text(face = "bold"),
          strip.text = element_text(size = 12, face = "bold")) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.5, alpha = 1)))
}
