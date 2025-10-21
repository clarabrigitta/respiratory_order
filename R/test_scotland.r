# load data
data <- read_csv("data/respiratory_scot_20250917.csv") %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>% 
  arrange(Pathogen, WeekBeginning) %>%
  group_by(Pathogen) %>%
  # mutate(epi_onset = rollapply(NumberCasesPerWeek, width = 5, FUN = function(v) all(diff(v) > 0), align = "right", fill = NA)) %>% # monotonic increase definition
  mutate(upweek = rollapply(NumberCasesPerWeek, width = 5, FUN = function(v) sum(diff(v) > 0) >= 3, align = "right", fill = NA)) %>% # 3/4 past weeks definition
  mutate(epi_onset = {
    x <- replace_na(upweek, FALSE)
    r <- rle(x)
    idx <- cumsum(r$lengths)
    onset_positions <- idx[r$values & r$lengths >= 4] - r$lengths[r$values & r$lengths >= 4] + 1
    seq_along(x) %in% onset_positions
  }) %>% # only keeping the first week if 4 consecutive upweeks to mark onset week
  mutate(forward = rollapply(NumberCasesPerWeek, width = 5, FUN = function(v) sum(diff(v) > 0) >= 3, align = "left", fill = NA)) %>%
  mutate(forward_onset = {
    x <- replace_na(forward, FALSE)
    r <- rle(x)
    idx <- cumsum(r$lengths)
    onset_positions <- idx[r$values & r$lengths >= 4] - r$lengths[r$values & r$lengths >= 4] + 1
    seq_along(x) %in% onset_positions
  }) %>% 
  mutate(growth = coalesce(
    (lead(NumberCasesPerWeek, 1) + lead(NumberCasesPerWeek, 2) + lead(NumberCasesPerWeek, 3)) >=
      3 * (lag(NumberCasesPerWeek, 1) + lag(NumberCasesPerWeek, 2) + lag(NumberCasesPerWeek, 3)),
    FALSE
  )) %>%
  mutate(growth_onset = {
    x <- replace_na(growth, FALSE)
    r <- rle(x)
    idx <- cumsum(r$lengths)
    onset_positions <- idx[r$values & r$lengths >= 4] - r$lengths[r$values & r$lengths >= 4] + 1
    seq_along(x) %in% onset_positions
  }) %>% 
  mutate(across(c(epi_onset, forward_onset, growth_onset), ~ ifelse(NumberCasesPerWeek <= 5 & .x == TRUE, FALSE, .x))) %>% 
  # mutate(epi_onset = NumberCasesPerWeek > lag(NumberCasesPerWeek, 4)) %>% # Zhao et al definition
  ungroup() %>% 
  filter(!Pathogen %in% c("Influenza - Type A (not subtyped)",
                          "Influenza - Type A or B",
                          "Influenza - Type A(H1N1)pdm09" ,
                          "Influenza - Type A(H3)",
                          "Mycoplasma pneumoniae" ))

# load data by age
data <- read_csv("data/cases_all_respiratory_pathogens_by_agegroup_sex_20251008.csv") %>% 
  filter(Sex == "Total", !AgeGroup %in% c("Total", "Unknown"), !Pathogen %in% c("Influenza (All)", "COVID-19")) %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) 

c("<1", "1 to 4", "5 to 14", "15 to 44", "45 to 64", "65 to 74", "75+")
ggplot() +
  geom_line(data = data %>% filter(AgeGroup %in% c("<1", "1 to 4", "5 to 14"),
                                   ISOyear >= 2018 & ISOyear <= 2022), aes(x = WeekBeginning, y = NumberCasesPerWeek, colour = Pathogen)) +
  scale_color_viridis(discrete = T, option = "D") +
  theme_bw() +
  facet_grid(~AgeGroup)

subdata <- data %>% filter(AgeGroup %in% c("<1", "1 to 4", "5 to 14"),
                               ISOyear >= 2018 & ISOyear <= 2022)

  subplot(lapply(split(subdata, subdata$AgeGroup), function(df) {
    plot_ly(
      data = df,
      x = ~WeekBeginning,
      y = ~NumberCasesPerWeek,
      legendgroup = ~Pathogen,
      color = ~Pathogen,
      colors = viridis(length(unique(df$Pathogen)), option = "D"),
      type = 'scatter',
      mode = 'lines',
      showlegend = ifelse(unique(df$AgeGroup) == "<1", TRUE, FALSE)) %>%
      layout(annotations = list(text = paste0(unique(df$AgeGroup)),
                                x = 0.5, y = 1, xref = 'paper', yref = 'paper',
                                showarrow = FALSE))}),
    nrows = 1, shareY = TRUE, shareX = TRUE)

# plot using ggplot
ggplot(data) +
  geom_line(aes(x = WeekBeginning, y = NumberCasesPerWeek, colour = Pathogen)) +
  theme_bw() +
  annotate("rect",
           xmin = as.Date("2020-03-23"), xmax = as.Date("2021-03-16"),
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "grey")

# plot using plotly
plot_ly() %>%
  add_trace(data = data, x = ~WeekBeginning, y = ~NumberCasesPerWeek, color = ~Pathogen,
            type = 'scatter', 
            mode = 'lines') %>% 
  add_trace(data = subset(data, epi_onset), x = ~WeekBeginning, y = ~NumberCasesPerWeek, color = ~Pathogen,
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 7, symbol = 'circle', line = list(width = 1, color = NA)),
            showlegend = FALSE) %>% 
  add_trace(data = subset(data, forward_onset), x = ~WeekBeginning, y = ~NumberCasesPerWeek, color = ~Pathogen,
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 7, symbol = 'x', line = list(width = 1, color = NA)),
            showlegend = FALSE) %>% 
  add_trace(data = subset(data, growth_onset), x = ~WeekBeginning, y = ~NumberCasesPerWeek, color = ~Pathogen,
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 7, symbol = 'square', line = list(width = 1, color = NA)),
            showlegend = FALSE) %>% 
  add_trace(data = contacts_daily, x = ~date, y = ~mean_contacts,
            yaxis = "y2",
            type = 'scatter', 
            mode = 'lines',
            color = 'contacts',
            line = list(color='black')) %>% 
  layout(shapes = list(list(type = "rect",
                            x0 = as.Date("2020-03-23"),
                            x1 = as.Date("2021-03-16"),
                            y0 = 0, y1 = 1,
                            xref = "x",
                            yref = "paper",  # spans the whole y-axis
                            fillcolor = "grey",
                            opacity = 0.2,
                            line = list(width = 0))), template = "plotly_white", 
         yaxis2 = list(overlaying = "y",
                       side = "right", range = c(0, 30)))
