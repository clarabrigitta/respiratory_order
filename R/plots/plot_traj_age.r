# Age-specific case trajectory panels, one saved figure per virus.
# Saves inst/plots/traj_<virus>_<date>.png and returns a named list of patchwork objects.
#
# results_traj : list of posterior trajectories, one element per virus (process_outputs.r)
# data         : weekly case counts, wide by age group (see plot_summary.r)
# combined     : AGE-STRATIFIED contact data with `date`, `agegp` and `mean_contacts`
#                (explore_contacts.r line 295)

plot_traj_age <- function(results_traj, data, combined, age_groups, pathogen_map,
                          outdir = here("inst", "plots")) {

  plots <- list()

  for (virus_name in names(results_traj)) {
    pathogen_name <- pathogen_map[[virus_name]]
    plots[[virus_name]] <- list()

    traj    <- results_traj[[virus_name]]
    subdata <- data %>%
      filter(Pathogen == pathogen_name) %>%
      select(WeekBeginning, all_of(age_groups)) %>%
      pivot_longer(-WeekBeginning, names_to = "age_group", values_to = "count") %>%
      mutate(age_group = factor(age_group, levels = age_groups))

    # HDI per age group across posterior samples
    week_dates <- as.Date(rownames(traj[[1]]))

    for (agegp in age_groups) {
      agedata <- do.call(rbind, lapply(traj, function(m) m[, agegp]))
      traj_hdi_age <- agedata %>%
        hdi() %>%
        rbind(mean = colMeans(agedata)) %>%
        t() %>%
        as.data.frame() %>%
        mutate(date = week_dates)

      subdata_age <- subdata %>% filter(age_group == agegp)

      local({
        agegp_local  <- agegp
        contact_max  <- max(combined$mean_contacts, na.rm = TRUE)
        # panel height must cover every series drawn, else geoms get clipped at the top
        y_max        <- max(c(subdata_age$count, traj_hdi_age$upper), na.rm = TRUE)
        scale_factor <- y_max / contact_max
        subcombined  <- combined %>% mutate(contacts_scaled = mean_contacts * scale_factor) %>% filter(agegp == agegp_local)

        plots[[virus_name]][[agegp_local]] <<- ggplot() +
          geom_ribbon(data = traj_hdi_age, aes(x = date, ymin = lower, ymax = upper, fill = "Model estimate"), alpha = 0.2) +
          geom_line(data = traj_hdi_age, aes(x = date, y = mean, colour = "Model estimate"), linewidth = 1) +
          geom_point(data = subdata_age, aes(x = WeekBeginning, y = count, colour = "Data"), size = 2) +
          geom_line(data = subcombined, aes(x = date, y = contacts_scaled), colour = "blue", linewidth = 0.5) +
          scale_color_manual(name = NULL, values = c("Data" = "black", "Model estimate" = "red")) +
          scale_fill_manual(name = NULL, values = c("Data" = "black", "Model estimate" = "red")) +
          scale_x_date(date_breaks = "3 month") +
          scale_y_continuous(name = "Count",
                             sec.axis = sec_axis(~ . / scale_factor, name = "Mean number of contacts",
                                                 breaks = seq(0, 20, 5))) +
          coord_cartesian(ylim = c(0, y_max)) +
          labs(x = "Weeks", title = agegp_local) +
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
  }

  figs <- list()

  for (virus_name in names(results_traj)) {

    fig_traj_age <- wrap_plots(plots[[virus_name]], ncol = 3) +
      plot_annotation(title = virus_name,
                      theme = theme(plot.title = element_text(size = 16, face = "bold"))) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")

    ggsave(filename = file.path(outdir, paste0("traj_", virus_name, "_", format(Sys.Date(), "%d%m%Y"), ".png")),
           plot = fig_traj_age, width = 15, height = 9, dpi = 300)

    figs[[virus_name]] <- fig_traj_age
  }

  figs
}
