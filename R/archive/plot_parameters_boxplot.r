# Panel plot of the range of each parameter estimate, one boxplot panel per parameter.
# Saves inst/plots/parameters_boxplot_<date>.png and returns the patchwork object.
#
# results : list of BayesianTools fits, one element per virus
# start   : first iteration kept (burn-in), as in process_parameters.r. Row 1 of each
#           chain is the sampler's start value, drawn from the prior.

plot_parameters_boxplot <- function(results, params, start = 3, thin = "auto",
                                    outdir = here("inst", "plots")) {

  data <- bind_rows(lapply(names(results), function(v) {
    getSample(results[[v]], start = start, thin = thin, parametersOnly = TRUE) %>%
      as.data.frame() %>%
      mutate(virus = v)})) %>%
    pivot_longer(-virus, names_to = "parameter", values_to = "value")

  plots <- lapply(params, function(p) {
    data %>%
      filter(parameter == p) %>%
      ggplot() +
      geom_boxplot(aes(x = virus, y = value, fill = virus),
                   outliers = FALSE, size = 0.3) +
      scale_fill_viridis_d(option = "H") +
      labs(x = "Virus", y = "Value", title = p, fill = "Virus") +
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))})

  fig_paramsbox <- wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "none")

  ggsave(filename = file.path(outdir, paste0("parameters_boxplot_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_paramsbox, width = 15, height = 9, dpi = 300)

  fig_paramsbox
}
