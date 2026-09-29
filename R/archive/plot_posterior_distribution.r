# Panel plot of the posterior density of each parameter estimate, coloured by virus.
# Saves inst/plots/posterior_distribution_<date>.png and returns the ggplot object.
#
# results : list of BayesianTools fits, one element per virus
# start   : first iteration kept (burn-in), as in process_parameters.r. Row 1 of each
#           chain is the sampler's start value, drawn from the prior; leaving it in
#           stretches the x-axis over a range the posterior never visits.
# thin    : "auto" keeps ~5000 draws per chain, which is ample for a density
# x_scale : "log10" puts every parameter on a log10 axis (imports are already stored
#           that way). Detection rates differ ~1000x between viruses, so on a linear
#           axis the tightly estimated ones collapse to a spike at zero.

plot_posterior_distribution <- function(results, start = 3, thin = "auto",
                                        x_scale = c("log10", "linear"),
                                        outdir = here("inst", "plots")) {

  x_scale <- match.arg(x_scale)

  data <- bind_rows(lapply(names(results), function(v) {
    getSample(results[[v]], start = start, thin = thin, parametersOnly = TRUE) %>%
      as.data.frame() %>%
      mutate(virus = v)})) %>%
    pivot_longer(-virus, names_to = "parameter", values_to = "value")

  if (x_scale == "log10") {
    data <- data %>%
      # log10_import is log10(imports) already, so it only needs renaming
      mutate(value = if_else(parameter == "log10_import", value,
                             log10(if_else(value > 0, value, NA_real_))),
             parameter = if_else(parameter == "log10_import", "import", parameter)) %>%
      filter(!is.na(value))
  }

  fig_postdist <- ggplot(data, aes(x = value, colour = virus, fill = virus)) +
    geom_density(alpha = 0.3, lwd = 0.8) +
    scale_colour_viridis_d() +
    scale_fill_viridis_d() +
    # few breaks: panels span very different magnitudes under scales = "free"
    scale_x_continuous(n.breaks = 5,
                       labels = scales::label_number(drop0trailing = TRUE)) +
    labs(x = if (x_scale == "log10") expression(log[10]~"(value)") else "Value",
         y = "Density") +
    facet_wrap(~ parameter, scales = "free", nrow = 2) +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          strip.text   = element_text(size = 12, face = "bold"))

  if (x_scale == "linear")
    fig_postdist <- fig_postdist +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(filename = file.path(outdir, paste0("posterior_distribution_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_postdist, width = 21, height = 6, dpi = 300)

  fig_postdist
}
