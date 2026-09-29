# MCMC traceplot panel: one row per parameter, one column per virus.
# Base graphics, so this writes straight to inst/plots/chains_<date>.png
# and returns the file path invisibly.
#
# results : list of BayesianTools fits, one element per virus

plot_chains <- function(results, params, n_chains = 4, n_subchains = 3,
                        outdir = here("inst", "plots")) {

  colours_vec <- viridis(n_chains * n_subchains)
  virus_names <- names(results)

  filename <- file.path(outdir, paste0("chains_", format(Sys.Date(), "%d%m%Y"), ".png"))

  png(filename = filename, width = 20, height = 16, units = "in", res = 300)

  on.exit(dev.off(), add = TRUE)

  par(mfrow = c(length(params), length(virus_names)),
      mar   = c(2, 5, 2, 1),
      oma   = c(0, 0, 3, 0))

  for (i in seq_along(params)) {
    for (j in seq_along(virus_names)) {
      chain_list <- as.mcmc.list(
        unlist(
          lapply(seq_len(n_chains), function(s)
            lapply(seq_len(n_subchains), function(k)
              as.mcmc(window(results[[virus_names[j]]][[s]][["chain"]][[k]], start = 10, drop = FALSE)[, params[i]])
            )
          ),
          recursive = FALSE
        )
      )
      coda::traceplot(chain_list, col = colours_vec, main = "")
      if (i == 1) mtext(virus_names[j], side = 3, line = 0.5, font = 2, cex = 0.8)
      if (j == 1) mtext(params[i],      side = 2, line = 3,   font = 2, cex = 0.8)
    }
  }

  invisible(filename)
}
