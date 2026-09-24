# Pairwise correlation plot for one virus, laid out like BayesianTools'
# correlationPlot() (default density = "smooth"):
#   diagonal : histogram of each parameter, scaled so the tallest bar is 2/3 of the panel
#   lower    : binned 2D density (IDPmisc colour ramp) with a lowess line, as ipanel.smooth
#   upper    : Pearson correlation, text size proportional to |r|
# Chains are pooled. Returns the ggplot object.
#
#   plot_correlation(results, "RSV")
#
# results : list of BayesianTools fits, one element per virus
# start   : first iteration kept (burn-in), as in process_parameters.r
# thin    : "auto" keeps ~5000 draws per chain, as correlationPlot() does
# bins    : bins per axis for the 2D density panels

plot_correlation <- function(results, virus, start = 3, thin = "auto",
                             method = "pearson", bins = 50) {

  stopifnot(virus %in% names(results))

  draws <- getSample(results[[virus]], start = start, thin = thin, parametersOnly = TRUE)
  par_names <- colnames(draws)
  rng <- unname(apply(draws, 2, range))
  lvl <- function(x) factor(x, levels = par_names)

  pairs_idx <- expand.grid(i = seq_along(par_names), j = seq_along(par_names))

  # lower panels (row below column): 2D density and lowess line
  lower <- pairs_idx[pairs_idx$i > pairs_idx$j, ]

  density_2d <- bind_rows(lapply(seq_len(nrow(lower)), function(k) {
    x <- draws[, lower$j[k]]; y <- draws[, lower$i[k]]
    bx <- seq(rng[1, lower$j[k]], rng[2, lower$j[k]], length.out = bins + 1)
    by <- seq(rng[1, lower$i[k]], rng[2, lower$i[k]], length.out = bins + 1)
    counts <- as.data.frame(table(x = cut(x, bx, include.lowest = TRUE, labels = FALSE),
                                  y = cut(y, by, include.lowest = TRUE, labels = FALSE)))
    counts %>%
      filter(Freq > 0) %>%
      # colour scaled within each panel, as ipanel.smooth does
      mutate(Freq = Freq / max(Freq),
             x = (bx[-1] + bx[-(bins + 1)])[as.integer(as.character(x))] / 2,
             y = (by[-1] + by[-(bins + 1)])[as.integer(as.character(y))] / 2,
             width = diff(bx)[1], height = diff(by)[1],
             col_par = par_names[lower$j[k]], row_par = par_names[lower$i[k]])}))

  smooth <- bind_rows(lapply(seq_len(nrow(lower)), function(k) {
    as.data.frame(lowess(draws[, lower$j[k]], draws[, lower$i[k]], f = 2/3, iter = 3)) %>%
      mutate(col_par = par_names[lower$j[k]], row_par = par_names[lower$i[k]])}))

  # diagonal: histogram, heights mapped into the row's y range (BT uses usr y = 0-1.5)
  hists <- bind_rows(lapply(seq_along(par_names), function(i) {
    h <- hist(draws[, i], plot = FALSE)
    lo <- rng[1, i]; span <- rng[2, i] - rng[1, i]
    data.frame(xmin = head(h$breaks, -1), xmax = h$breaks[-1],
               ymin = lo, ymax = lo + (h$counts / max(h$counts)) / 1.5 * span,
               col_par = par_names[i], row_par = par_names[i])}))

  # upper panels (row above column): correlation text in the middle of the panel
  upper <- pairs_idx[pairs_idx$i < pairs_idx$j, ]
  cors <- cor(draws, method = method)
  cor_text <- data.frame(
    x = colMeans(rng)[upper$j], y = colMeans(rng)[upper$i],
    r = cors[cbind(upper$i, upper$j)],
    col_par = par_names[upper$j], row_par = par_names[upper$i])

  # blank points at each panel's corners, so free scales cover the full range
  # of both parameters even in text-only and histogram panels
  frame <- bind_rows(lapply(seq_len(nrow(pairs_idx)), function(k) {
    data.frame(x = rng[, pairs_idx$j[k]], y = rng[, pairs_idx$i[k]],
               col_par = par_names[pairs_idx$j[k]], row_par = par_names[pairs_idx$i[k]])}))

  facet_levels <- function(d) mutate(d, col_par = lvl(col_par), row_par = lvl(row_par))

  ggplot() +
    geom_blank(data = facet_levels(frame), aes(x = x, y = y)) +
    geom_tile(data = facet_levels(density_2d),
              aes(x = x, y = y, width = width, height = height, fill = Freq)) +
    geom_path(data = facet_levels(smooth), aes(x = x, y = y), colour = "black", linewidth = 0.4) +
    geom_rect(data = facet_levels(hists),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "blue4", colour = "white", linewidth = 0.1) +
    geom_text(data = facet_levels(cor_text),
              aes(x = x, y = y, label = sprintf("%.2f", r), size = abs(r))) +
    scale_fill_gradientn(colours = IDPmisc::IDPcolorRamp(100), limits = c(0, 1), guide = "none") +
    # text size linear in |r| like correlationPlot(), floored so r near 0 stays legible
    scale_radius(range = c(2.5, 8), limits = c(0, 1), guide = "none") +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3) +
    facet_grid(row_par ~ col_par, scales = "free") +
    labs(title = virus, x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3),
          panel.spacing = unit(0.15, "lines"),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(size = 8, face = "bold"),
          strip.text.y = element_text(angle = 0))
}
