# Panel plot of total (all ages summed) case trajectories for all pathogens.
# Saves inst/plots/traj_all_<date>.png and returns the patchwork object.
#
# results_traj : list of posterior trajectories, one element per virus (process_outputs.r)
# data         : weekly case counts, wide by age group (see plot_summary.r)
# combined     : NON-age-stratified contact data with `date` and `mean_contacts`
#                (explore_contacts.r line 271)

plot_traj_total <- function(results_traj, data, combined, age_groups, pathogen_map,
                            outdir = here("inst", "plots")) {

  plots <- list()

  for (virus_name in names(results_traj)) {
    pathogen_name <- pathogen_map[[virus_name]]

    traj <- results_traj[[virus_name]]
    subdata <- data %>%
      filter(Pathogen == pathogen_name) %>%
      select(WeekBeginning, all_of(age_groups)) %>%
      mutate(total = rowSums(across(all_of(age_groups)))) %>%
      select(-all_of(age_groups))

    traj_all <- do.call(cbind, lapply(traj, rowSums)) %>% as.data.frame() %>% t()
    traj_hdi <- traj_all %>%
      hdi() %>%
      rbind(mean = colMeans(traj_all)) %>%
      t() %>%
      as.data.frame() %>%
      mutate(date = unique(subdata$WeekBeginning))

    local({
      scale_factor <- max(subdata$total, na.rm = TRUE) / max(combined$mean_contacts, na.rm = TRUE)
      subcombined <- combined %>% mutate(contacts_scaled = mean_contacts * scale_factor)
      contact_breaks <- seq(0, ceiling(max(combined$mean_contacts)), by = 2)

      plots[[virus_name]] <<- ggplot() +
        geom_ribbon(data = traj_hdi, aes(x = date, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
        geom_line(data = traj_hdi, aes(x = date, y = mean, color = "Model estimate"), linewidth = 1) +
        geom_line(data = subdata, aes(x = WeekBeginning, y = total, color = "Data"), alpha = 0.7, linewidth = 1) +
        geom_line(data = subcombined, aes(x = date, y = contacts_scaled), colour = "blue", lty = 2, linewidth = 1) +
        scale_color_manual(name = NULL, values = c("Data" = "black", "Model estimate" = "red")) +
        scale_x_date(date_breaks = "3 month") +
        scale_y_continuous(name = "Count",
                           sec.axis = sec_axis(~ . / scale_factor, name = "Mean number of contacts",
                                               breaks = contact_breaks)) +
        labs(x = "Weeks", title = virus_name) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 11),
              plot.title = element_text(size = 12, face = "bold"),
              axis.title.y.right = element_text(colour = "blue"),
              axis.text.y.right = element_text(colour = "blue"),
              axis.line.y.right = element_line(colour = "blue"),
              axis.ticks.y.right = element_line(colour = "blue"))
    })
  }

  fig_traj <- wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(filename = file.path(outdir, paste0("traj_all_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_traj, width = 15, height = 9, dpi = 300)

  fig_traj
}
