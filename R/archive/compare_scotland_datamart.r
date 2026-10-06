library(readxl)
library(purrr)
library(viridis)

# combine all sheets into one dataframe
sheet <- excel_sheets("inst/data/LSHTM_data_request.xlsx")

datamart <- map(set_names(sheet), function(nm) {
  read_excel("inst/data/LSHTM_data_request.xlsx", sheet = nm) %>%
    rename_with(~ "P_count", contains("P_count")) %>%
    rename_with(~ "N_count", contains("N_count")) 
})

datamart <- bind_rows(datamart, .id = "pathogen") %>% 
  mutate(week_no = as.Date(paste0(week_no, "1"), format = "%Y%W%u"),
         total_count = P_count + N_count,
         positivity = P_count/(N_count+P_count)*100,
         check = ifelse(P_count>N_count, TRUE, FALSE)) %>% 
  mutate(age_group = factor(age_group, levels = c("0-4", "5-14", "15-44", "45-64", "65+")))

# count time series
datamart_age <- ggplot(datamart 
       %>% filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))
       ) +
  geom_line(aes(x = week_no, y = P_count, colour = pathogen)) +
  # geom_vline(xintercept = as.Date("2020-03-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2020-08-11"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-01-05"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-04-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-08-09"), linetype = "dashed") + 
  theme_bw() +
  scale_colour_manual(values = pathogen_cols, na.value = "grey70") +
  # scale_x_date(date_breaks = "1 month") +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Number of Cases per Week", colour = "Pathogen", title = "Datamart") +
  facet_wrap(~age_group, nrow = 1, scales = "free_y")

ggsave(filename = here("inst", "plots", "datamart", paste0("datamart_count_age_freey", ".png")),
       plot = datamart_age, width = 12, height = 5, dpi = 300)

datamart_pathogen <- ggplot(datamart 
       %>% filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))
) +
  geom_line(aes(x = week_no, y = P_count, colour = age_group)) +
  # geom_vline(xintercept = as.Date("2020-03-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2020-08-11"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-01-05"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-04-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-08-09"), linetype = "dashed") + 
  theme_bw() +
  scale_colour_viridis(option = "H", discrete = T, begin = 0.1, end = 0.85) +
  # scale_x_date(date_breaks = "1 month") +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Number of Cases per Week", colour = "Age Group", title = "Datamart") +
  facet_wrap(~pathogen, nrow = 2, scales = "free_y")

ggsave(filename = here("inst", "plots", "datamart", paste0("datamart_count_pathogen_freey", ".png")),
       plot = datamart_pathogen, width = 16, height = 10, dpi = 300)

# stacked bar plot (bottom bar is positive counts, top is greyed out number of tests)
datamart %>%
  filter(week_no >= as.Date("2021-03-23") & week_no <= as.Date("2022-03-02")) %>%
  pivot_longer(cols = c(P_count, N_count), names_to = "type", values_to = "count") %>%
  ggplot() +
  geom_col(aes(x = week_no, y = count, fill = type)) +
  theme_bw() +
  scale_fill_manual(values = c("N_count" = "lightgrey", "P_count" = "lightgreen")) +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Week Beginning", y = "Count", fill = NULL) +
  facet_grid(pathogen~age_group)

# positivity time series
pos_age <- ggplot(datamart 
       %>% filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))
) +
  geom_line(aes(x = week_no, y = positivity, colour = pathogen)) +
  # geom_vline(xintercept = as.Date("2020-03-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2020-08-11"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-01-05"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-04-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-08-09"), linetype = "dashed") + 
  theme_bw() +
  scale_colour_manual(values = pathogen_cols, na.value = "grey70") +
  # scale_x_date(date_breaks = "1 month") +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Positivity (%)", colour = "Pathogen") +
  facet_grid(~age_group)

ggsave(filename = here("inst", "plots", "datamart", paste0("datamart_positivity_age", ".png")),
       plot = pos_age, width = 12, height = 5, dpi = 300)

pos_pathogen <- ggplot(datamart %>% 
                         filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))) +
  geom_line(aes(x = week_no, y = positivity, colour = age_group)) +
  # geom_vline(xintercept = as.Date("2020-03-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2020-08-11"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-01-05"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-04-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-08-09"), linetype = "dashed") + 
  theme_bw() +
  scale_colour_viridis(option = "H", discrete = T, begin = 0.1, end = 0.85) +
  # scale_x_date(date_breaks = "1 month") +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Positivity (%)", colour = "Age Group") +
  facet_grid(~pathogen)

ggsave(filename = here("inst", "plots", "datamart", paste0("datamart_positivity_pathogen", ".png")),
       plot = pos_pathogen, width = 16, height = 5, dpi = 300)

# testing effort against positivity time series
plot_df <- datamart %>%
  filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))

coeff <- max(plot_df$total_count, na.rm = TRUE) / max(plot_df$positivity, na.rm = TRUE)

test_pos <- ggplot(plot_df) +
  geom_col(aes(x = week_no, y = total_count, fill = pathogen), alpha = 0.5, fill = "lightgrey") +
  geom_line(aes(x = week_no, y = positivity * coeff, colour = pathogen)) +
  theme_bw() +
  scale_colour_manual(values = pathogen_cols, na.value = "grey70", guide = "none") +
  scale_y_continuous(
    name = "Number of Tests per Week",
    sec.axis = sec_axis(~ . / coeff, name = "Positivity (%)")
  ) +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(x = "Week Beginning", fill = "Pathogen") +
  facet_grid(pathogen ~ age_group)

ggsave(filename = here("inst", "plots", "datamart", paste0("datamart_tests_positivity", ".png")),
       plot = test_pos, width = 12, height = 10, dpi = 300)

# compare scottish and datamart data
data <- read_csv("inst/data/cases_all_respiratory_pathogens_by_agegroup_sex_20251008.csv") %>%
  filter(Sex == "Total",
         !AgeGroup %in% c("Total", "Unknown"),
         !Pathogen %in% c("Influenza (All)", "COVID-19", "Mycoplasma pneumoniae")) %>%
  mutate(WeekBeginning = as.Date(as.character(WeekBeginning), format = "%Y%m%d")) %>%
  filter(WeekBeginning >= fit_start,
         WeekBeginning <= fit_end) %>%
  pivot_wider(names_from = AgeGroup, values_from = NumberCasesPerWeek) %>%
  mutate(`0 to 4` = `<1` + `1 to 4`,
         `65+` = `65 to 74` + `75+`) %>%
  select(-c(`<1`, `1 to 4`, `65 to 74`, `75+`)) %>%
  arrange(WeekBeginning) %>% 
  mutate(Pathogen = case_when(Pathogen == "Adenovirus" ~ "Adeno",
                              Pathogen == "Influenza A" ~ "FluA",
                              Pathogen == "Influenza B" ~ "FluB",
                              Pathogen == "Parainfluenza (Any Type)" ~ "Parainfluenza",
                              Pathogen == "Rhinovirus" ~ "Rhino_Entero",
                              Pathogen == "Seasonal coronavirus" ~ "Seasonal_Corona",
                              .default = Pathogen)) %>% 
  pivot_longer(cols = c(`5 to 14`:`65+`), names_to = "age_group", values_to = "count") %>% 
  mutate(age_group = case_when(age_group == "0 to 4" ~ "0-4",
                               age_group == "5 to 14" ~ "5-14",
                               age_group == "15 to 44" ~ "15-44",
                               age_group == "45 to 64" ~ "45-64",
                              .default = age_group)) %>% 
  mutate(age_group = factor(age_group, levels = c("0-4", "5-14", "15-44", "45-64", "65+"))) %>% 
  rename(pathogen = Pathogen, week_no = WeekBeginning, P_count = count) %>% 
  mutate(source = "Scotland") %>% 
  select(source, week_no, pathogen, age_group, P_count) %>% 
  bind_rows(datamart %>% mutate(source = "Datamart") %>% select(source, week_no, pathogen, age_group, P_count))

scot_age <- ggplot(data %>% filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))) +
  geom_line(aes(x = week_no, y = P_count, colour = pathogen)) +
  # geom_vline(xintercept = as.Date("2020-03-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2020-08-11"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-01-05"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-04-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-08-09"), linetype = "dashed") + 
  theme_bw() +
  scale_colour_manual(values = pathogen_cols, na.value = "grey70") +
  # scale_x_date(date_breaks = "1 month") +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Number of Cases per Week", colour = "Pathogen", title = "Scottish") +
  facet_wrap(~age_group, nrow = 1, scales = "free_y")

ggsave(filename = here("inst", "plots", "datamart", paste0("scotland_count_age_freey", ".png")),
       plot = scot_age, width = 12, height = 5, dpi = 300)

scot_pathogen <- ggplot(data
       %>% filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))
) +
  geom_line(aes(x = week_no, y = P_count, colour = age_group)) +
  # geom_vline(xintercept = as.Date("2020-03-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2020-08-11"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-01-05"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-04-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-08-09"), linetype = "dashed") + 
  theme_bw() +
  scale_colour_viridis(option = "H", discrete = T, begin = 0.1, end = 0.85) +
  # scale_x_date(date_breaks = "1 month") +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Number of Cases per Week", colour = "Age Group", title = "Scottish") +
  facet_wrap(~pathogen, nrow = 2, scales = "free_y")

ggsave(filename = here("inst", "plots", "datamart", paste0("scotland_count_pathogen_freey", ".png")),
       plot = scot_pathogen, width = 16, height = 10, dpi = 300)

both <- ggplot(data %>% filter(week_no >= as.Date("2020-03-23") & week_no <= as.Date("2022-03-02"))) +
  geom_line(aes(x = week_no, y = P_count, colour = pathogen, linetype = source)) +
  # geom_vline(xintercept = as.Date("2020-03-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2020-08-11"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-01-05"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-04-26"), linetype = "dashed") + 
  # geom_vline(xintercept = as.Date("2021-08-09"), linetype = "dashed") + 
  theme_bw() +
  scale_colour_manual(values = pathogen_cols, na.value = "grey70") +
  # scale_x_date(date_breaks = "1 month") +
  scale_x_date(date_breaks = "6 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Number of Cases per Week", colour = "Pathogen") +
  facet_grid(pathogen~age_group)

ggsave(filename = here("inst", "plots", "datamart", paste0("datamart_scotland_count", ".png")),
       plot = both, width = 12, height = 11, dpi = 300)
