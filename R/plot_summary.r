results <- readRDS(file = here("inst", "outdata", "parameters_rcpp_sus_det_imm_13072026"))
results_traj <- readRDS(file = here("inst", "outdata", "traj_rcpp_sus_det_imm_13072026"))

# panel plot of trajectories for all pathogens ----
data <- read_csv("inst/data/cases_all_respiratory_pathogens_by_agegroup_sex_20251008.csv") %>% 
  filter(Sex == "Total", !AgeGroup %in% c("Total", "Unknown"), !Pathogen %in% c("Influenza (All)", "COVID-19", "Mycoplasma pneumoniae")) %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>% 
  filter(WeekBeginning >= as.Date("2020-03-23"),
         WeekBeginning <= as.Date("2022-03-02")) %>% 
  pivot_wider(names_from = AgeGroup, values_from = NumberCasesPerWeek) %>% 
  mutate(`0 to 4` = `<1` + `1 to 4`,
         `65+` = `65 to 74` + `75+`) %>% 
  select(-c(`<1`, `1 to 4`, `65 to 74`, `75+`)) %>% 
  arrange(WeekBeginning)

# load age group bands
age_groups  <- c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")
# age_colours <- setNames(viridis(5), age_groups)
# need to also load combined from explore_contacts.r either total or age-stratified (lines 271 or 295)

# total ----
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

ggsave(filename = here("inst", "plots", paste0("traj_all_", format(Sys.Date(), "%d%m%Y"), ".png")), plot = fig_traj, width = 15, height = 9, dpi = 300)

# age-specific ----
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
      # colour_val   <- age_colours[[agegp_local]]
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

for (virus_name in names(results_traj)) {
  
  fig_traj_age <- wrap_plots(plots[[virus_name]], ncol = 3) +
    plot_annotation(title = virus_name,
                    theme = theme(plot.title = element_text(size = 16, face = "bold"))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(filename = here("inst", "plots", paste0("traj_", virus_name, "_", format(Sys.Date(), "%d%m%Y"), ".png")), plot = fig_traj_age, width = 15, height = 9, dpi = 300)
  
}


# MCMC traceplot panel for all viruses and parameters ----
n_chains    <- 4
n_subchains <- 3
colours_vec <- viridis(n_chains * n_subchains)
params      <- c("detection_0to4", "detection_5to14", "detection_15to44", "detection_45to64", "detection_65plus", "sus_dist", "imm_duration")
virus_names <- names(results)

png(
  filename = here("inst", "plots", paste0("chains_", format(Sys.Date(), "%d%m%Y"), ".png")),
  width = 20, height = 12, units = "in", res = 300
)

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

dev.off()

# panel plot of posterior distribution of parameter estimates ----
data <- bind_rows(lapply(names(results), function(v) {
  getSample(results[[v]], skip = 10) %>% 
    as.data.frame() %>% 
    mutate(virus = v)})) %>% 
  pivot_longer(-virus, names_to = "parameter", values_to = "value")

fig_postdist <- ggplot(data, aes(x = value, colour = virus, fill = virus)) +
  geom_density(alpha = 0.3, lwd = 0.8) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  # scale_x_continuous(limits = c(100, 500)) +
  labs(x = "Value", y = "Density") +
  facet_wrap(~ parameter, scales = "free", nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        strip.text   = element_text(size = 12, face = "bold"))

ggsave(filename = here("inst", "plots", paste0("posterior_distribution_", format(Sys.Date(), "%d%m%Y"), ".png")), plot = fig_postdist, width = 18, height = 6, dpi = 300)

# panel plot of range of parameter estimates ----
data <- bind_rows(lapply(names(results), function(v) {
  getSample(results[[v]], skip = 10) %>% 
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

ggsave(filename = here("inst", "plots", paste0("parameters_boxplot_", format(Sys.Date(), "%d%m%Y"), ".png")), plot = fig_paramsbox, width = 15, height = 9, dpi = 300)

# plot Rt trajectories ----
plots <- list()

for (virus_name in names(results_rt)) {
  pathogen_name <- pathogen_map[[virus_name]]
  
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

fig_traj <- wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = here("inst", "plots", paste0("rt_", format(Sys.Date(), "%d%m%Y"), ".png")), plot = fig_traj, width = 15, height = 9, dpi = 300)

# plot susceptible S(t) by age group for each virus, with mean contacts on 2nd axis ----

# map the 9 model age groups onto the 5 reporting bands (same as incidence aggregation)
s_age_map <- list(`0 to 4`   = 1,
                  `5 to 14`  = 2:3,
                  `15 to 44` = 4:6,
                  `45 to 64` = 7:8,
                  `65+`      = 9)
day_dates <- as.Date("2020-03-23") + times  # daily date for each ODE time step

# per-age-group mean contacts (shared across all virus panels)
contacts_age <- combined %>%
  transmute(date, age_group = factor(agegp, levels = age_groups), mean_contacts) %>%
  drop_na(age_group)

plots <- list()

for (virus_name in names(results)) {
  n <- which(vapply(combinations, function(x) x$name, character(1)) == virus_name)
  
  posterior <- getSample(results[[virus_name]], start = 2, thin = 100)
  
  # accumulate S (aggregated to 5 age bands) across posterior samples, then average
  s_sum <- matrix(0, nrow = length(times), ncol = length(s_age_map))
  
  for (r in seq_len(nrow(posterior))) {
    sus_dist     <- posterior[r, 6]
    imm_duration <- posterior[r, 7]
    
    R_init <- floor(tots * (1 - sus_dist))
    y0_mat <- rbind(S = tots - 5 - R_init, E = rep(0, 9), I = rep(5, 9), R = R_init)
    y0 <- setNames(as.vector(y0_mat),
                   paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4)))
    
    ode_out <- seirs_rcpp(
      y0    = y0,
      times = times,
      parms = list(
        sigma = 1 / combinations[[n]]$inc_period,
        gamma = 1 / combinations[[n]]$inf_period,
        omega = 1 / imm_duration,
        p_inf = combinations[[n]]$p_inf,
        mu_b  = daily_births,
        mu_d  = hazard_death
      ),
      contacts_prepped = contacts_prepped
    )
    
    S_all <- ode_out[, seq(2, 34, by = 4)]  # S1..S9 (S state column for each age group)
    s_sum <- s_sum + vapply(s_age_map,
                            function(cols) rowSums(S_all[, cols, drop = FALSE]),
                            numeric(length(times)))
  }
  
  s_mean <- s_sum / nrow(posterior)
  colnames(s_mean) <- names(s_age_map)
  
  s_df <- as.data.frame(s_mean) %>%
    mutate(date = day_dates) %>%
    pivot_longer(-date, names_to = "age_group", values_to = "S") %>%
    mutate(age_group = factor(age_group, levels = age_groups))
  
  local({
    scale_factor <- max(s_df$S, na.rm = TRUE) / max(contacts_age$mean_contacts, na.rm = TRUE)
    sub_contacts <- contacts_age %>% mutate(contacts_scaled = mean_contacts * scale_factor)
    contact_breaks <- seq(0, ceiling(max(contacts_age$mean_contacts, na.rm = TRUE)), by = 2)
    
    plots[[virus_name]] <<- ggplot() +
      geom_line(data = s_df, aes(x = date, y = S, colour = age_group),
                linewidth = 0.8) +
      geom_line(data = sub_contacts, aes(x = date, y = contacts_scaled, colour = age_group),
                linetype = "dashed", linewidth = 0.5) +
      scale_colour_viridis_d(option = "D", end = 0.9) +
      scale_x_date(date_breaks = "3 month") +
      scale_y_continuous(name = "Susceptible (S)",
                         sec.axis = sec_axis(~ . / scale_factor,
                                             name = "Mean number of contacts",
                                             breaks = contact_breaks)) +
      labs(x = "Date", colour = "Age group", title = virus_name) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text  = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12, face = "bold"))
  })
}

fig_susceptible <- wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = here("inst", "plots", paste0("susceptible_contacts_", format(Sys.Date(), "%d%m%Y"), ".png")),
       plot = fig_susceptible, width = 15, height = 9, dpi = 300)

# plot TOTAL susceptible S(t) for each virus, with mean contacts on 2nd axis ----
# Requires model_rcpp.r sourced (seirs_rcpp, combinations, times) and
# contacts_prepped, tots, daily_births, hazard_death available (see fit_model_rcpp.r).
# Requires the NON-age-stratified `combined` (with `mean_contacts`) from explore_contacts.r line 271.

day_dates <- as.Date("2020-03-23") + times  # daily date for each ODE time step

plots <- list()

for (virus_name in names(results)) {
  n <- which(vapply(combinations, function(x) x$name, character(1)) == virus_name)
  
  posterior <- getSample(results[[virus_name]], start = 2, thin = 100)
  
  # accumulate total S (all age groups summed) across posterior samples, then average
  s_sum <- numeric(length(times))
  
  for (r in seq_len(nrow(posterior))) {
    sus_dist     <- posterior[r, 6]
    imm_duration <- posterior[r, 7]
    
    R_init <- floor(tots * (1 - sus_dist))
    y0_mat <- rbind(S = tots - 5 - R_init, E = rep(0, 9), I = rep(5, 9), R = R_init)
    y0 <- setNames(as.vector(y0_mat),
                   paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4)))
    
    ode_out <- seirs_rcpp(
      y0    = y0,
      times = times,
      parms = list(
        sigma = 1 / combinations[[n]]$inc_period,
        gamma = 1 / combinations[[n]]$inf_period,
        omega = 1 / imm_duration,
        p_inf = combinations[[n]]$p_inf,
        mu_b  = daily_births,
        mu_d  = hazard_death
      ),
      contacts_prepped = contacts_prepped
    )
    
    S_all <- ode_out[, seq(2, 34, by = 4)]  # S1..S9
    s_sum <- s_sum + rowSums(S_all)          # total susceptible across all age groups
  }
  
  s_df <- data.frame(date = day_dates, S = s_sum / nrow(posterior))
  
  local({
    scale_factor <- max(s_df$S, na.rm = TRUE) / max(combined$mean_contacts, na.rm = TRUE)
    sub_contacts <- combined %>% mutate(contacts_scaled = mean_contacts * scale_factor)
    contact_breaks <- seq(0, ceiling(max(combined$mean_contacts, na.rm = TRUE)), by = 2)
    
    plots[[virus_name]] <<- ggplot() +
      geom_line(data = s_df, aes(x = date, y = S), colour = "darkgreen", linewidth = 0.8) +
      geom_line(data = sub_contacts, aes(x = date, y = contacts_scaled),
                colour = "blue", linetype = "dashed", linewidth = 0.6) +
      scale_x_date(date_breaks = "3 month") +
      scale_y_continuous(name = "Total susceptible (S)",
                         sec.axis = sec_axis(~ . / scale_factor,
                                             name = "Mean number of contacts",
                                             breaks = contact_breaks)) +
      labs(x = "Date", title = virus_name) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text  = element_text(size = 10),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12, face = "bold"),
            axis.title.y.right = element_text(colour = "blue"),
            axis.text.y.right  = element_text(colour = "blue"),
            axis.line.y.right  = element_line(colour = "blue"),
            axis.ticks.y.right = element_line(colour = "blue"))
  })
}

fig_susceptible_total <- wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = here("inst", "plots", paste0("susceptible_total_", format(Sys.Date(), "%d%m%Y"), ".png")),
       plot = fig_susceptible_total, width = 15, height = 9, dpi = 300)

