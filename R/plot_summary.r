# panel plot of trajectories for all pathogens ----
plots <- list()

for (virus_name in names(results_traj)) {
  pathogen_name <- pathogen_map[[virus_name]]
  
  traj <- results_traj[[virus_name]]
  subdata <- data %>% 
    filter(Pathogen == pathogen_name) 
  
  traj_all <- do.call(cbind, traj) %>%  as.data.frame() %>%  t()
  traj_hdi <- traj_all %>% 
    hdi() %>% 
    rbind(mean = colMeans(traj_all)) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(date = unique(subdata$WeekBeginning))
  
  local({
    scale_factor <- max(subdata$NumberCasesPerWeek, na.rm = TRUE) / max(combined$mean_contacts, na.rm = TRUE)
    subcombined <- combined %>% mutate(contacts_scaled = mean_contacts * scale_factor)
    contact_breaks <- seq(0, ceiling(max(combined$mean_contacts)), by = 2)

    plots[[virus_name]] <<- ggplot() +
      geom_ribbon(data = traj_hdi, aes(x = date, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
      geom_line(data = traj_hdi, aes(x = date, y = mean, color = "Model estimate"), linewidth = 1) +
      geom_line(data = subdata, aes(x = WeekBeginning, y = NumberCasesPerWeek, color = "Data"), alpha = 0.7, linewidth = 1) +
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

# MCMC traceplot panel for all viruses and parameters ----
colours_vec <- viridis(3)
params <- c("sus_dist", "detection_rate", "imm_duration")
virus_names <- names(results)

par(mfrow = c(length(params), length(virus_names)),
    mar   = c(2, 5, 2, 1),
    oma   = c(0, 0, 3, 0))

for (i in seq_along(params)) {
  for (j in seq_along(virus_names)) {
    traceplot(results[[virus_names[j]]][["chain"]][, params[i]],
              col = vec_colours,main = "")
    if (i == 1) mtext(virus_names[j], side = 3, line = 0.5, font = 2, cex = 0.8)
    if (j == 1) mtext(params[i],      side = 2, line = 3,   font = 2, cex = 0.8)
  }
}

dev.off() # to avoid any downstream issues with plotting after par()

# panel plot of posterior distribution of parameter estimates ----
data <- bind_rows(lapply(names(results), function(v) {
  getSample(results[[v]]) %>% 
    as.data.frame() %>% 
    mutate(virus = v)})) %>% 
  pivot_longer(-virus, names_to = "parameter", values_to = "value")

fig_postdist <- ggplot(data, aes(x = value, colour = virus, fill = virus)) +
  geom_density(alpha = 0.3, lwd = 0.8) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  labs(x = "Value", y = "Density") +
  facet_wrap(~ parameter, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        strip.text   = element_text(size = 12, face = "bold"))

# panel plot of range of parameter estimates ----
data <- bind_rows(lapply(names(results), function(v) {
  getSample(results[[v]]) %>% 
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

# plot Rt trajectories ----
plots <- list()

for (virus_name in names(results_traj)) {
  pathogen_name <- pathogen_map[[virus_name]]
  
  rt_traj <- results_rt[[virus_name]]
  
  rt_hdi <- do.call(cbind, lapply(rt_traj, `[[`, "rt")) %>%
    t() %>% 
    hdi() %>%
    rbind(mean = colMeans(do.call(cbind, lapply(rt_traj, `[[`, "rt")) %>% t())) %>%
    t() %>% 
    as.data.frame() %>%
    mutate(date = rt_trajectories[[1]]$date)
  
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

fig_traj <- wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
