# Panel plot of Rt trajectories for all viruses.
# Saves inst/plots/rt_<date>.png and returns the patchwork object.
#
# results_rt : list of Rt trajectories, one element per virus (process_outputs.r)

plot_rt <- function(results_rt, outdir = here("inst", "plots")) {

  plots <- list()

  for (virus_name in names(results_rt)) {

    rt_traj <- results_rt[[virus_name]]

    rt_hdi <- do.call(cbind, lapply(rt_traj, `[[`, "rt")) %>%
      t() %>%
      hdi() %>%
      rbind(mean = colMeans(do.call(cbind, lapply(rt_traj, `[[`, "rt")) %>% t())) %>%
      t() %>%
      as.data.frame() %>%
      mutate(date = rt_traj[[1]]$date)

    local({

      plots[[virus_name]] <<- ggplot() +
        geom_hline(yintercept = 1, lty = 2, colour = "red") +
        geom_ribbon(data = rt_hdi, aes(x = date, ymin = lower, ymax = upper), fill = "black", alpha = 0.2) +
        geom_line(data = rt_hdi, aes(x = date, y = mean), colour = "black", linewidth = 1) +
        scale_x_date(date_breaks = "3 month") +
        scale_y_continuous(name = "R_t") +
        labs(x = "Date", y = "R_t", title = virus_name) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              plot.title = element_text(size = 12, face = "bold"))
    })
  }

  fig_rt <- wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(filename = file.path(outdir, paste0("rt_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_rt, width = 15, height = 9, dpi = 300)

  fig_rt
}
