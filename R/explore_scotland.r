# load Scottish data by pathogen, marking epidemic onset ----
data <- read_csv("inst/data/respiratory_scot_20250917.csv") %>% 
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

# plot Scottish data by pathogen ----
## ggplot
ggplot(data %>% filter(WeekBeginning >= as.Date("2020-03-23") & WeekBeginning <= as.Date("2022-03-02"))) +
  geom_line(aes(x = WeekBeginning, y = NumberCasesPerWeek, colour = Pathogen)) +
  theme_bw() +
  scale_colour_viridis_d(option = "H", 
                         begin = 0,
                         end = 1) +
  scale_x_date(date_breaks = "1 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      legend.text=element_text(size=12),
      legend.title=element_text(size=14))
  # xlim(as.Date("2020-03-23"), as.Date("2021-03-16")) +
  # annotate("rect",
  #          # xmin = as.Date("2020-03-23"), xmax = as.Date("2021-03-16"),
  #          ymin = -Inf, ymax = Inf,
  #          alpha = 0.2, fill = "grey")

## plotly
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
  # add_trace(data = contacts_daily, x = ~date, y = ~mean_contacts,
  #           yaxis = "y2",
  #           type = 'scatter', 
  #           mode = 'lines',
  #           color = 'contacts',
  #           line = list(color='black')) %>% 
  # add_trace(data = out, x = ~date, y = ~inc,
  #           yaxis = "y2",
  #           type = 'scatter',
  #           mode = 'lines',
  #           color = 'model',
  #           line = list(color='black')) %>%
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
                       side = "right", range = c(0, 100000)))

# load Scottish data by pathogen and age ----
data <- read_csv("inst/data/cases_all_respiratory_pathogens_by_agegroup_sex_20251008.csv") %>% 
  filter(Sex == "Total", !AgeGroup %in% c("Total", "Unknown"), !Pathogen %in% c("Influenza (All)", "COVID-19")) %>% 
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) 

# c("<1", "1 to 4", "5 to 14", "15 to 44", "45 to 64", "65 to 74", "75+") # vector of age ranges

# plot Scottish data by pathogen and age ----
## ggplot
ggplot() +
  geom_line(data = data %>% filter(AgeGroup %in% c("<1", "1 to 4", "5 to 14"),
                                   ISOyear >= 2020 & ISOyear <= 2021), aes(x = WeekBeginning, y = NumberCasesPerWeek, colour = Pathogen)) +
  scale_color_viridis(discrete = T, option = "D") +
  theme_bw() +
  facet_grid(~AgeGroup)

## plotly
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
  
  
# load Scottish mid-population estimates ----
## data source: https://www.nrscotland.gov.uk/publications/population-estimates-time-series-data/
scot_population <- read_excel("inst/data/mid-year-population-estimates-time-series-data.xlsx", sheet = "Table 1", skip = 5) %>% 
    filter(`Area name` == "Scotland", Sex == "Persons", Year %in% c(2020:2022)) %>% 
    select(-`All Ages`) %>% 
    pivot_longer(cols = 5:95, names_to = "age", values_to = "n") %>% 
    mutate(agegp = case_when(age %in% c("0", "1", "2", "3", "4") ~ "0-4",
                             age %in% c("5", "6", "7", "8", "9", "10", "11") ~ "5-11",
                             age %in% c("12", "13", "14", "15", "16", "17") ~ "12-17",
                             age %in% c("18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29") ~ "18-29",
                             age %in% c("30", "31", "32", "33", "34", "35", "36", "37", "38", "39") ~ "30-39",
                             age %in% c("40", "41", "42", "43", "44", "45", "46", "47", "48", "49") ~ "40-49",
                             age %in% c("50", "51", "52", "53", "54", "55", "56", "57", "58", "59") ~ "50-59",
                             age %in% c("60", "61", "62", "63", "64", "65", "66", "67", "68", "69") ~ "60-69",
                             age %in% c(as.character(c(70:89)), "90 and over") ~"70+")) %>% 
    group_by(agegp, Year) %>% 
    summarise(sum_n = sum(n)) %>% 
    ungroup() %>% 
    group_by(agegp) %>% 
    summarise(average_n = floor(mean(sum_n))) %>% # taking the average between 2020, 2021 and 2022
    ungroup() %>% 
    mutate(agegp = factor(agegp, levels = c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+")))
    
    