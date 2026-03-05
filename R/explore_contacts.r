# CoMix data (all ages, up to 2021) ----
comix <- read_csv("inst/data/contact_matrices_9_periods.csv") %>% 
  group_by(Period, `Participant age`) %>%
  summarise(contacts = sum(`mean contacts`, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(Period) %>% 
  summarise(mean_contacts = mean(contacts, na.rm = TRUE)) %>% 
  rename(period = Period) %>% 
  mutate(period = 1:9,
         period = factor(period, 
                         levels = c(1:9), 
                         labels = c("Lockdown 1", "Lockdown 1 easing", "Relaxed restrictions", "School reopening", "Lockdown 2", "Lockdown 2 easing", "Christmas", "Lockdown 3", "Lockdown 3 + schools")))

contacts_daily <- data.frame(date = seq(from = as.Date("23-03-2020", format = "%d-%m-%Y"), to   = as.Date("16-03-2021", format = "%d-%m-%Y"), by   = "day")) %>% 
  mutate(time = as.numeric(1:nrow(.))) %>% 
  mutate(period = case_when(
    date >= as.Date("23-03-2020", format = "%d-%m-%Y") & date <= as.Date("03-06-2020", format = "%d-%m-%Y") ~ 1,
    date >= as.Date("04-06-2020", format = "%d-%m-%Y") & date <= as.Date("29-07-2020", format = "%d-%m-%Y") ~ 2,
    date >= as.Date("30-07-2020", format = "%d-%m-%Y") & date <= as.Date("03-09-2020", format = "%d-%m-%Y") ~ 3,
    date >= as.Date("04-09-2020", format = "%d-%m-%Y") & date <= as.Date("26-10-2020", format = "%d-%m-%Y") ~ 4,
    date >= as.Date("27-10-2020", format = "%d-%m-%Y") & date <= as.Date("02-12-2020", format = "%d-%m-%Y") ~ 5, # actually 05-11-2025
    date >= as.Date("03-12-2020", format = "%d-%m-%Y") & date <= as.Date("19-12-2020", format = "%d-%m-%Y") ~ 6,
    date >= as.Date("20-12-2020", format = "%d-%m-%Y") & date <= as.Date("02-01-2021", format = "%d-%m-%Y") ~ 7,
    date >= as.Date("03-01-2021", format = "%d-%m-%Y") & date <= as.Date("08-03-2021", format = "%d-%m-%Y") ~ 8, # actually 05-01-2025
    date >= as.Date("09-03-2021", format = "%d-%m-%Y") & date <= as.Date("16-03-2021", format = "%d-%m-%Y") ~ 9)) %>%
  mutate(period = factor(period, 
                         levels = c(1:9), 
                         labels = c("Lockdown 1", "Lockdown 1 easing", "Reduced restrictions", "Schools open", "Lockdown 2", "Lockdown 2 easing", "Christmas", "Lockdown 3", "Lockdown 3 with schools open"))) %>% 
  left_join(comix)

ggplot() +
  geom_line(data = contacts_daily, aes(x = date, y = mean_contacts)) +
  geom_line(data = data %>% filter(!Pathogen %in% c("Rhinovirus", "Adenovirus")), aes(x = WeekBeginning, y = NumberCasesPerWeek, colour = Pathogen)) +
  theme_bw()

# CoMix data (age-stratified, up to 2021) ----
comix_kids <- read_csv("inst/data/contact_matrices_9_periods.csv") %>% 
  group_by(Period, `Participant age`) %>%
  summarise(contacts = sum(`mean contacts`, na.rm = TRUE)) %>% 
  # filter(`Participant age` %in% c("0-4", "5-11", "12-17")) %>% 
  mutate(`Participant age` = factor(`Participant age`, levels = c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+"))) %>% 
  rename(period = Period) %>% 
  mutate(period = recode(period, 
                         "1. Lockdown 1" = "Lockdown 1", 
                         "2. Lockdown 1 easing" = "Lockdown 1 easing", 
                         "3. Relaxed restrictions" = "Relaxed restrictions", 
                         "4. School reopening" = "School reopening", 
                         "5. Lockdown 2" = "Lockdown 2", 
                         "6. Lockdown 2 easing" = "Lockdown 2 easing", 
                         "7. Christmas" = "Christmas", 
                         "8. Lockdown 3" = "Lockdown 3", 
                         "9. Lockdown 3 + schools" = "Lockdown 3 + schools"))


contacts_daily_kids <- data.frame(date = seq(from = as.Date("23-03-2020", format = "%d-%m-%Y"), to   = as.Date("16-03-2021", format = "%d-%m-%Y"), by   = "day")) %>% 
  mutate(time = as.numeric(1:nrow(.))) %>% 
  mutate(period = case_when(
    date >= as.Date("23-03-2020", format = "%d-%m-%Y") & date <= as.Date("03-06-2020", format = "%d-%m-%Y") ~ 1,
    date >= as.Date("04-06-2020", format = "%d-%m-%Y") & date <= as.Date("29-07-2020", format = "%d-%m-%Y") ~ 2,
    date >= as.Date("30-07-2020", format = "%d-%m-%Y") & date <= as.Date("03-09-2020", format = "%d-%m-%Y") ~ 3,
    date >= as.Date("04-09-2020", format = "%d-%m-%Y") & date <= as.Date("26-10-2020", format = "%d-%m-%Y") ~ 4,
    date >= as.Date("27-10-2020", format = "%d-%m-%Y") & date <= as.Date("02-12-2020", format = "%d-%m-%Y") ~ 5, # actually 05-11-2025
    date >= as.Date("03-12-2020", format = "%d-%m-%Y") & date <= as.Date("19-12-2020", format = "%d-%m-%Y") ~ 6,
    date >= as.Date("20-12-2020", format = "%d-%m-%Y") & date <= as.Date("02-01-2021", format = "%d-%m-%Y") ~ 7,
    date >= as.Date("03-01-2021", format = "%d-%m-%Y") & date <= as.Date("08-03-2021", format = "%d-%m-%Y") ~ 8, # actually 05-01-2025
    date >= as.Date("09-03-2021", format = "%d-%m-%Y") & date <= as.Date("16-03-2021", format = "%d-%m-%Y") ~ 9)) %>%
  mutate(period = factor(period, 
                         levels = c(1:9), 
                         labels = c("Lockdown 1", "Lockdown 1 easing", "Relaxed restrictions", "School reopening", "Lockdown 2", "Lockdown 2 easing", "Christmas", "Lockdown 3", "Lockdown 3 + schools"))) %>% 
  left_join(comix_kids)

ggplot() +
  geom_line(data = contacts_daily_kids, aes(x = date, y = contacts, colour = `Participant age`)) +
  scale_color_viridis(discrete = T, option = "D") +
  theme_bw()

# CoMix data (all ages, up to 2022) ----

survey <- read_csv("inst/data/CoMix_uk_sday.csv") %>% 
  select(-X)

participants <- read_csv("inst/data/CoMix_uk_participant_common.csv") %>% 
  select(part_id, part_age)

contacts <- read_csv("inst/data/CoMix_uk_contact_common.csv") %>% 
  group_by(part_id) %>% 
  summarise(count = n())

date_period <- data.frame(date = seq(from = as.Date("23-03-2020", format = "%d-%m-%Y"), to   = as.Date("02-03-2022", format = "%d-%m-%Y"), by   = "day")) %>% 
  mutate(time = as.numeric(1:nrow(.)),
         mmyyyy = format(date, "%m/%Y")) %>% 
  mutate(period = case_when(
    date >= as.Date("23-03-2020", format = "%d-%m-%Y") & date <= as.Date("03-06-2020", format = "%d-%m-%Y") ~ 1,
    date >= as.Date("04-06-2020", format = "%d-%m-%Y") & date <= as.Date("29-07-2020", format = "%d-%m-%Y") ~ 2,
    date >= as.Date("30-07-2020", format = "%d-%m-%Y") & date <= as.Date("03-09-2020", format = "%d-%m-%Y") ~ 3,
    date >= as.Date("04-09-2020", format = "%d-%m-%Y") & date <= as.Date("26-10-2020", format = "%d-%m-%Y") ~ 4,
    date >= as.Date("27-10-2020", format = "%d-%m-%Y") & date <= as.Date("02-12-2020", format = "%d-%m-%Y") ~ 5, # actually 05-11-2025
    date >= as.Date("03-12-2020", format = "%d-%m-%Y") & date <= as.Date("19-12-2020", format = "%d-%m-%Y") ~ 6,
    date >= as.Date("20-12-2020", format = "%d-%m-%Y") & date <= as.Date("02-01-2021", format = "%d-%m-%Y") ~ 7,
    date >= as.Date("03-01-2021", format = "%d-%m-%Y") & date <= as.Date("08-03-2021", format = "%d-%m-%Y") ~ 8, # actually 05-01-2025
    date >= as.Date("09-03-2021", format = "%d-%m-%Y") & date <= as.Date("16-03-2021", format = "%d-%m-%Y") ~ 9,
    TRUE ~ 10)) %>% 
  mutate(period = factor(period,
                         levels = c(1:10), 
                         labels = c("Lockdown 1", "Lockdown 1 easing", "Relaxed restrictions", "School reopening", "Lockdown 2", "Lockdown 2 easing", "Christmas", "Lockdown 3", "Lockdown 3 + schools", "Post lockdown")))

combined <- survey %>% 
  left_join(participants, join_by(part_id)) %>% 
  left_join(contacts, join_by(part_id)) %>% 
  mutate(count = replace_na(count, 0),
         sday_id = as.Date(sday_id, format = "%Y.%m.%d"),
         # mmyyyy = format(sday_id, "%m/%Y"),
         mmyyyy = as.yearmon(sday_id),
         period = case_when(
           sday_id >= as.Date("23-03-2020", format = "%d-%m-%Y") & sday_id <= as.Date("03-06-2020", format = "%d-%m-%Y") ~ 1,
           sday_id >= as.Date("04-06-2020", format = "%d-%m-%Y") & sday_id <= as.Date("29-07-2020", format = "%d-%m-%Y") ~ 2,
           sday_id >= as.Date("30-07-2020", format = "%d-%m-%Y") & sday_id <= as.Date("03-09-2020", format = "%d-%m-%Y") ~ 3,
           sday_id >= as.Date("04-09-2020", format = "%d-%m-%Y") & sday_id <= as.Date("26-10-2020", format = "%d-%m-%Y") ~ 4,
           sday_id >= as.Date("27-10-2020", format = "%d-%m-%Y") & sday_id <= as.Date("02-12-2020", format = "%d-%m-%Y") ~ 5, # actually 05-11-2025
           sday_id >= as.Date("03-12-2020", format = "%d-%m-%Y") & sday_id <= as.Date("19-12-2020", format = "%d-%m-%Y") ~ 6,
           sday_id >= as.Date("20-12-2020", format = "%d-%m-%Y") & sday_id <= as.Date("02-01-2021", format = "%d-%m-%Y") ~ 7,
           sday_id >= as.Date("03-01-2021", format = "%d-%m-%Y") & sday_id <= as.Date("08-03-2021", format = "%d-%m-%Y") ~ 8, # actually 05-01-2025
           sday_id >= as.Date("09-03-2021", format = "%d-%m-%Y") & sday_id <= as.Date("16-03-2021", format = "%d-%m-%Y") ~ 9,
           TRUE ~ 10)) %>% 
  # filter(sday_id < as.Date("31-03-2021", format = "%d-%m-%Y")) %>% 
  group_by(mmyyyy) %>% 
  summarise (mean_contacts = mean(count, na.rm = TRUE)) %>% 
  mutate(scaled = (mean_contacts - min(mean_contacts))/(max(mean_contacts) - min(mean_contacts))) %>%
  # mutate(period = factor(period, 
  #                        levels = c(1:10), 
  #                        labels = c("Lockdown 1", "Lockdown 1 easing", "Relaxed restrictions", "School reopening", "Lockdown 2", "Lockdown 2 easing", "Christmas", "Lockdown 3", "Lockdown 3 + schools", "Post lockdown"))) %>% 
  # mutate(part_age = factor(part_age,
  #                          levels = c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70-120"))) %>% 
  # drop_na(part_age) %>% 
  left_join(date_period %>% mutate(mmyyyy = as.yearmon(date)), join_by(mmyyyy))

ggplot() +
  geom_line(data = combined, aes(x = date, y = scaled)) +
  scale_color_viridis(discrete = T, option = "D") +
  theme_bw()
  
# Age-stuctured contact matrices (using CoMix and socialmixr) ----

## load comix data from zenodo
comix_uk <- get_survey("https://doi.org/10.5281/zenodo.4905745")

## load participant file to fill missing participant age
participants <- read_csv("inst/data/CoMix_uk_participant_common.csv") %>% 
  select(part_id, part_age) %>% 
  mutate(part_age = recode(part_age, "Under 1" = "0-1")) %>% 
  separate(part_age, into = c("part_age_est_min", "part_age_est_max"), sep = "-", convert = TRUE)

## fill in missing participant age
comix_uk$participants <- comix_uk$participants %>%
  left_join(participants, by = "part_id")

## save complete UK comix data
# saveRDS(comix_uk, "inst/data/comix_uk.rds")

# contact_matrix <- contact_matrix(comix_uk, age.limits = c(0, 5, 12, 18, 30, 40, 50, 60, 70)) # test for entire survey period
contact_matrix <- contact_matrix(comix_uk, age.limits = c(0, 5, 11, 15, 25, 35, 45, 55, 65)) # test for entire survey period

## define fortnight variable in participant data
participants <- comix_uk$participants %>%
  mutate(date = as.Date(sday_id, format = "%Y.%m.%d"),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2)))) %>% 
  group_by(fortnight) %>% 
  mutate(start_date = min(date),
         end_date = max (date),
         mid_date = start_date + floor((end_date - start_date)/2)) %>% 
  ungroup() %>% 
  group_split(mid_date) # separate data by fortnights

## define fortnight variable in contact data
survey <- read_csv("inst/data/CoMix_uk_sday.csv") %>% 
  select(part_id, sday_id)

contacts <- comix_uk$contacts %>%
  left_join(survey, join_by(part_id)) %>% 
  mutate(date = as.Date(sday_id, format = "%Y.%m.%d"),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2)))) %>% 
  group_by(fortnight) %>% 
  mutate(start_date = min(date),
         end_date = max (date),
         mid_date = start_date + floor((end_date - start_date)/2)) %>% 
  ungroup() %>% 
  group_split(mid_date) # separate data by fortnights


## create surveys split by fortnight
fortnight_survey <- list()

for (i in 1:length(participants)) {
  fortnight_survey[[i]] <- survey(participants[[i]], contacts[[i]])
}

## create contact matrix for each fortnight survey
fortnight_matrix <- list()

for (i in 1:length(participants)) {
  fortnight_matrix[[i]] <- contact_matrix(fortnight_survey[[i]],
                                          age.limits = c(0, 5, 11, 15, 25, 35, 45, 55, 65))
                                          # age.limits = c(0, 1, 5, 10, 15, 25, 35, 45, 55, 65))
                                          # age.limits = c(0, 5, 12, 18, 30, 40, 50, 60, 70))
}

## identify which matrices contain missing data
na_fortnight <- logical(length(fortnight_matrix))

for (i in seq_along(fortnight_matrix)) {  
  m <- fortnight_matrix[[i]]$matrix
  
  if (anyNA(m)) {
    na_fortnight[i] <- TRUE
  }
}

## example plot of single fortnightly matrix
fortnight_matrix[[20]]$matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("part_age_group") %>% 
  rename("0-4" = "[0,5)", "5-10" = "[5,11)", "11-14" = "[11,15)", "15-24" = "[15,25)", "25-34" = "[25,35)", "35-44" = "[35,45)", "45-55" = "[45,55)", "55-64" = "[55,65)") %>% 
  mutate(part_age_group = case_match(part_age_group,
                                     "[0,5)" ~ "0-4",
                                     "[5,11)" ~ "5-10",
                                     "[11,15)" ~ "11-14",
                                     "[15,25)" ~ "15-24",
                                     "[25,35)" ~ "25-34",
                                     "[35,45)" ~ "35-44",
                                     "[45,55)" ~ "45-55",
                                     "[55,65)" ~ "55-64",
                                     "65+" ~ "65+")) %>%
  pivot_longer(-part_age_group, names_to = "contact_age_group", values_to = "mean_contacts") %>% 
  mutate(part_age_group = factor(part_age_group, levels = c("0-4", "5-10", "11-14", "15-24", "25-34", "35-44", "45-55", "55-64", "65+")),
         contact_age_group = factor(contact_age_group, levels = c("0-4", "5-10", "11-14", "15-24", "25-34", "35-44", "45-55", "55-64", "65+"))) %>% 
  ggplot(aes(x = part_age_group, y = contact_age_group, fill = mean_contacts)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D") +
  theme_bw() +
  labs(x = "Participant Age Group", y = "Contact Age Group", fill = "Mean Contacts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

# processing for first 3 fortnights and last fortnight containing missing participant data ----

fortnight_matrix[[1]][["matrix"]] <- fortnight_matrix[[4]][["matrix"]]
fortnight_matrix[[2]][["matrix"]] <- fortnight_matrix[[4]][["matrix"]]
fortnight_matrix[[3]][["matrix"]] <- fortnight_matrix[[4]][["matrix"]]
fortnight_matrix[[52]][["matrix"]] <- fortnight_matrix[[51]][["matrix"]]

# Fortnightly average number of contacts per age group ----
survey <- read_csv("inst/data/CoMix_uk_sday.csv") %>% 
  select(-X)

participants <- read_csv("inst/data/CoMix_uk_participant_common.csv") %>% 
  select(part_id, part_age)

contacts <- read_csv("inst/data/CoMix_uk_contact_common.csv") %>% 
  group_by(part_id) %>% 
  summarise(count = n())

dates <- data.frame(date = seq(from = as.Date("23-03-2020", format = "%d-%m-%Y"), to = as.Date("02-03-2022", format = "%d-%m-%Y"), by = "day")) %>% 
  mutate(fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date))

combined <- survey %>% 
  left_join(participants, join_by(part_id)) %>% 
  left_join(contacts, join_by(part_id)) %>% 
  mutate(count = replace_na(count, 0),
         date = as.Date(sday_id, format = "%Y.%m.%d"),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date)) %>% 
  group_by(quarter, part_age) %>%
  summarise(mean_contacts = mean(count, na.rm = TRUE)) %>% 
  mutate(part_age = factor(part_age,
                           levels = c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70-120"))) %>% 
  drop_na(part_age) %>% 
  left_join(dates, by = "quarter")

ggplot() +
  geom_line(data = combined, aes(x = date, y = mean_contacts, colour = part_age)) +
  scale_colour_viridis_d(option = "G", end = 0.8) + 
  scale_x_date(date_breaks = "1 month") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) + 
  labs(x = "Date", y = "Mean Number of Contacts", colour = "Age Group")

# POLYMOD data ----

## load POLYMOD data from zenodo
polymod <- get_survey("https://doi.org/10.5281/zenodo.3874557")

## process POLYMOD data
polymod[["participants"]] <- polymod[["participants"]] %>% filter(country == "United Kingdom") # filter for participants in UK
uk_part <- polymod[["participants"]] %>% select(part_id) %>% pull() # select UK participants id
polymod[["contacts"]] <- polymod[["contacts"]] %>% filter(part_id %in% uk_part, cnt_school == FALSE) # filter out school contacts

polymod_matrix <- contact_matrix(polymod, age.limits = c(0, 5, 12, 18, 30, 40, 50, 60, 70))

# imputation/processing for missing participant data ----

## compute dominant eigenvalue of 3 matrices missing data
missing1 <- fortnight_matrix[[1]][["matrix"]][4:9, 4:9]
dom1 <- max(Mod(eigen(missing1)$values))

missing2 <- fortnight_matrix[[2]][["matrix"]][4:9, 4:9]
dom2 <- max(Mod(eigen(missing2)$values))

missing3 <- fortnight_matrix[[3]][["matrix"]][4:9, 4:9]
dom3 <- max(Mod(eigen(missing3)$values))

## compute dominant eigenvalue of corresponding POLYMOD matrix
polymod_subset <- polymod_matrix[["matrix"]][4:9, 4:9]
dom_polymod <- max(Mod(eigen(polymod_subset)$values))

## compute scaling factor for 3 CoMix matrices
q1 <- dom1/dom_polymod
q2 <- dom2/dom_polymod
q3 <- dom3/dom_polymod

## replace missing CoMiX data with scaled POLYMOD
fortnight_matrix[[1]][["matrix"]][1:3, ] <- polymod_matrix[["matrix"]][1:3, ] * q1
fortnight_matrix[[1]][["matrix"]][4:9, 1:3] <- polymod_matrix[["matrix"]][4:9, 1:3] * q1

fortnight_matrix[[2]][["matrix"]][1:3, ] <- polymod_matrix[["matrix"]][1:3, ] * q2
fortnight_matrix[[2]][["matrix"]][4:9, 1:3] <- polymod_matrix[["matrix"]][4:9, 1:3] * q2

fortnight_matrix[[3]][["matrix"]][1:3, ] <- polymod_matrix[["matrix"]][1:3, ] * q3
fortnight_matrix[[3]][["matrix"]][4:9, 1:3] <- polymod_matrix[["matrix"]][4:9, 1:3] * q3
