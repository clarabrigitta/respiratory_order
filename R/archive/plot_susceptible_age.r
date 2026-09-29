# Panel plot of susceptible S(t) by age group for each virus, with mean contacts on 2nd axis.
# Saves inst/plots/susceptible_contacts_<date>.png and returns the patchwork object.
#
# Requires model_rcpp.r sourced (seirs_rcpp) and the model machinery from
# fit_model_rcpp.r passed in below.
# combined : AGE-STRATIFIED contact data with `date`, `agegp` and `mean_contacts`
#            (explore_contacts.r line 295)

plot_susceptible_age <- function(results, combined, combinations, age_groups,
                                 times, tots, daily_births, hazard_death, contacts_prepped,
                                 start_date = as.Date("2020-03-23"),
                                 outdir = here("inst", "plots")) {

  # map the 9 model age groups onto the 5 reporting bands (same as incidence aggregation)
  s_age_map <- list(`0 to 4`   = 1,
                    `5 to 14`  = 2:3,
                    `15 to 44` = 4:6,
                    `45 to 64` = 7:8,
                    `65+`      = 9)
  day_dates <- start_date + times  # daily date for each ODE time step

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

  ggsave(filename = file.path(outdir, paste0("susceptible_contacts_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_susceptible, width = 15, height = 9, dpi = 300)

  fig_susceptible
}
