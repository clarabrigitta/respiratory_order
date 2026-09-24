# ---------------------------------------------------------------------------
# Prior vs posterior, one facet per parameter, for a single virus.
#   plot_prior_posterior(results, combinations, "RSV")
# ---------------------------------------------------------------------------

# Draw from the prior exactly as createPrior() does in fit_model_rcpp.r:
#   params 1-10, 12 and 13 ~ Uniform(lb, ub); 13 is p_inf, with U(0, 1)
#   param 11 (immunity)    ~ Lognormal(log(imm_period), imm_sd), NOT truncated
#                            (the prior uses [0, Inf]; lb[11]/ub[11] are unused)
sample_prior <- function(combo, par_names, n = 1e5) {
  idx   <- list(det = 1:5, sus = 6:10, imm = 11, imp = 12, pinf = 13)
  u_idx <- c(idx$det, idx$sus, idx$imp, idx$pinf)
  
  out <- matrix(NA_real_, nrow = n, ncol = length(par_names),
                dimnames = list(NULL, par_names))
  for (j in u_idx) out[, j] <- runif(n, combo$lb[j], combo$ub[j])
  out[, idx$imm] <- rlnorm(n, meanlog = log(combo$imm_period), sdlog = combo$imm_sd)
  out
}

plot_prior_posterior <- function(results, combinations, virus,
                                 start = 3, thin = "auto", n_prior = 1e5,
                                 x_scale = c("log10", "linear"),
                                 normalise = TRUE) {
  
  x_scale <- match.arg(x_scale)
  stopifnot(virus %in% names(results))
  
  combo <- combinations[[which(vapply(combinations, `[[`, "", "name") == virus)]]
  if (length(combo) == 0) stop("no entry named '", virus, "' in combinations")
  
  post <- getSample(results[[virus]], start = start, thin = thin,
                    parametersOnly = TRUE)
  par_names <- colnames(post)
  
  long <- bind_rows(
    as.data.frame(post) %>% mutate(distribution = "posterior"),
    as.data.frame(sample_prior(combo, par_names, n_prior)) %>%
      mutate(distribution = "prior")
  ) %>%
    pivot_longer(-distribution, names_to = "parameter", values_to = "value") %>%
    mutate(parameter = factor(parameter, levels = par_names),
           distribution = factor(distribution, levels = c("prior", "posterior")))
  
  if (x_scale == "log10") {
    long <- long %>%
      # log10_import is already on a log10 scale, so it only needs relabelling
      mutate(value = if_else(parameter == "log10_import", value,
                             log10(if_else(value > 0, value, NA_real_)))) %>%
      filter(is.finite(value))
    levels(long$parameter)[levels(long$parameter) == "log10_import"] <- "import"
  }
  
  # one KDE per parameter x distribution, so each keeps its own bandwidth
  dens <- long %>%
    group_by(parameter, distribution) %>%
    reframe({
      y  <- value[!is.na(value)]
      bw <- 1.06 * min(sd(y), IQR(y) / 1.34) * length(y)^-0.2
      d  <- density(y, width = 4 * bw, n = 512)
      data.frame(x = d$x, y = d$y)
    }) %>%
    group_by(parameter, distribution) %>%
    mutate(y = if (normalise) y / max(y) else y) %>%
    ungroup()
  
  ggplot(dens, aes(x = x, y = y, colour = distribution, fill = distribution)) +
    geom_area(alpha = 0.25, position = "identity", linewidth = 0) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = c(prior = "grey55", posterior = "#2c7fb8")) +
    scale_fill_manual(values   = c(prior = "grey70", posterior = "#2c7fb8")) +
    scale_x_continuous(n.breaks = 5) +
    facet_wrap(~parameter, scales = "free", nrow = 3, ncol = 5) +
    labs(title = virus,
         x = if (x_scale == "log10") expression(log[10]~"(value)") else "Value",
         y = if (normalise) "Density (scaled to peak)" else "Density") +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
}
