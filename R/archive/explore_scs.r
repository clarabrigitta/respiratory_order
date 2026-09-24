# Scottish Contact Survey (SCS) contact matrices ----

dates <- data.frame(date = seq(from = as.Date("23-03-2020", format = "%d-%m-%Y"), to = as.Date("02-03-2022", format = "%d-%m-%Y"), by = "day")) %>% 
  mutate(fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date))

## age groups used in the survey
age_groups <- c("0-4", "5-12", "13-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+")

## load scs contact matrices and define fortnight variable
## note: survey waves are weekly up to 17-03-2022 and fortnightly after
scs <- read_csv("inst/data/scottish-contact-survey-contact-matrices.csv") %>%
  select(date = DateCode, part_age_group = `Participant Age Group`, cnt_age_group = `Contact Age Group`, mean_contacts = Value) %>%
  mutate(date = as.Date(date),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2)))) %>%
  mutate(part_age_group = case_match(part_age_group,
                                     "0-4 years" ~ "0-4",
                                     "5-12 years" ~ "5-12",
                                     "13-17 years" ~ "13-17",
                                     "18-29 years" ~ "18-29",
                                     "30-39 years" ~ "30-39",
                                     "40-49 years" ~ "40-49",
                                     "50-59 years" ~ "50-59",
                                     "60-69 years" ~ "60-69",
                                     "70 years and over" ~ "70+"),
         cnt_age_group = case_match(cnt_age_group,
                                    "0-4 years" ~ "0-4",
                                    "5-12 years" ~ "5-12",
                                    "13-17 years" ~ "13-17",
                                    "18-29 years" ~ "18-29",
                                    "30-39 years" ~ "30-39",
                                    "40-49 years" ~ "40-49",
                                    "50-59 years" ~ "50-59",
                                    "60-69 years" ~ "60-69",
                                    "70 years and over" ~ "70+"))


# Fortnightly average number of contacts per age group ----

## total contacts reported by each participant age group (row sums of each fortnightly matrix)
contacts_scs <- scs %>%
  mutate(part_age_group = factor(part_age_group, levels = age_groups)) %>% 
  group_by(date, fortnight, part_age_group) %>% 
  summarise(mean_contacts = sum(mean_contacts, na.rm = TRUE)) %>% # na.rm as child-child contacts are not reported
  ungroup() %>%
  group_by(fortnight, part_age_group) %>% 
  summarise(mean_contacts = mean(mean_contacts, na.rm = TRUE)) %>% 
  right_join(dates, by = "fortnight", relationship = "many-to-many") # keeps every day in the comix window, NA where scs has no wave

ggplot() +
  geom_line(data = contacts_scs, aes(x = date, y = mean_contacts, colour = part_age_group)) +
  scale_colour_viridis_d(option = "H", end = 0.8, na.translate = FALSE) + # days outside the scs period carry no age group
  scale_x_date(date_breaks = "2 months") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Date", y = "Mean Number of Contacts", colour = "Age Group")

# Fortnightly average number of contacts in total ----

## average across participant age groups
contacts_total <- contacts_agegp %>%
  group_by(fortnight, mid_date, date) %>%
  summarise(mean_contacts = mean(mean_contacts)) %>% # no na.rm, so days outside the scs period stay NA and break the line
  ungroup() %>%
  arrange(date)

ggplot() +
  geom_line(data = contacts_total, aes(x = date, y = mean_contacts), colour = "blue", lty = 2) +
  scale_x_date(date_breaks = "2 months") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Date", y = "Mean Number of Contacts")



## average survey waves within each fortnight (weekly waves only, one wave per fortnight from 2022)
scs_fortnight <- scs %>%
  group_by(fortnight, mid_date, part_age_group, cnt_age_group) %>%
  summarise(mean_contacts = mean(mean_contacts, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(nesting(fortnight, mid_date), part_age_group, cnt_age_group) %>% # child-child contacts (0-17 with 0-17) are not reported in the scs
  arrange(mid_date) %>%
  group_split(mid_date) # separate data by fortnights

## create contact matrix for each fortnight
fortnight_matrix <- list()

for (i in 1:length(scs_fortnight)) {
  fortnight_matrix[[i]] <- scs_fortnight[[i]] %>%
    select(part_age_group, cnt_age_group, mean_contacts) %>%
    pivot_wider(names_from = cnt_age_group, values_from = mean_contacts) %>%
    column_to_rownames("part_age_group") %>%
    as.matrix()
}

names(fortnight_matrix) <- dates$fortnight

## identify which matrices contain missing data
## note: all fortnights are missing the child-child block (0-17 with 0-17) and nothing else
na_fortnight <- logical(length(fortnight_matrix))

for (i in seq_along(fortnight_matrix)) {
  m <- fortnight_matrix[[i]]

  if (anyNA(m)) {
    na_fortnight[i] <- TRUE
  }
}

## example plot of single fortnightly matrix
fortnight_matrix[[20]] %>%
  as.data.frame() %>%
  rownames_to_column("part_age_group") %>%
  pivot_longer(-part_age_group, names_to = "cnt_age_group", values_to = "mean_contacts") %>%
  mutate(part_age_group = factor(part_age_group, levels = age_groups),
         cnt_age_group = factor(cnt_age_group, levels = age_groups)) %>%
  ggplot(aes(x = part_age_group, y = cnt_age_group, fill = mean_contacts)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D") +
  theme_bw() +
  labs(x = "Participant Age Group", y = "Contact Age Group", fill = "Mean Contacts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

## plot all fortnightly matrices
scs_fortnight %>%
  bind_rows() %>%
  ggplot(aes(x = part_age_group, y = cnt_age_group, fill = mean_contacts)) +
  geom_tile() +
  facet_wrap(~ mid_date) +
  scale_fill_viridis_c(option = "D") +
  theme_bw() +
  labs(x = "Participant Age Group", y = "Contact Age Group", fill = "Mean Contacts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
