# Panel plot of TOTAL susceptible S(t) for each virus, with mean contacts on 2nd axis.
# Saves inst/plots/susceptible_total_<date>.png and returns the patchwork object.
#
# Requires model_rcpp.r sourced (seirs_rcpp) and the model machinery from
# fit_model_rcpp.r passed in below.
# combined : NON-age-stratified contact data with `date` and `mean_contacts`
#            (explore_contacts.r line 271)

plot_susceptible_total <- function(results, combined, combinations,
                                   times, tots, daily_births, hazard_death, contacts_prepped,
                                   start_date = as.Date("2020-03-23"),
                                   outdir = here("inst", "plots")) {

  day_dates <- start_date + times  # daily date for each ODE time step

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

  ggsave(filename = file.path(outdir, paste0("susceptible_total_", format(Sys.Date(), "%d%m%Y"), ".png")),
         plot = fig_susceptible_total, width = 15, height = 9, dpi = 300)

  fig_susceptible_total
}
