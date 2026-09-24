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

## process population data for matrices
uk_pop <- read_excel("inst/data/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx", sheet = "Estimates", skip = 16) %>% 
  filter(`Region, subregion, country or area *` == "United Kingdom",
         Year %in% 2020:2022) %>% 
  select(Year, c(`0`:`100+`)) %>% 
  mutate(across(`0`:`100+`, as.numeric)) %>% 
  pivot_longer(cols = `0`:`100+`, names_to = "age", values_to = "pop") %>%
  pivot_wider(names_from = Year, values_from = pop) %>%
  mutate(mean = rowMeans(across(`2020`:`2022`), na.rm = TRUE))  

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
                                          age.limits = c(0, 5, 11, 15, 25, 35, 45, 55, 65),
                                          symmetric = TRUE)
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

# processing for first 3 fortnights and last fortnight containing missing participant data

fortnight_matrix[[1]][["matrix"]] <- fortnight_matrix[[4]][["matrix"]]
fortnight_matrix[[2]][["matrix"]] <- fortnight_matrix[[4]][["matrix"]]
fortnight_matrix[[3]][["matrix"]] <- fortnight_matrix[[4]][["matrix"]]
fortnight_matrix[[52]][["matrix"]] <- fortnight_matrix[[51]][["matrix"]]

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

## not age-stratified example
combined <- survey %>% 
  left_join(participants, join_by(part_id)) %>% 
  left_join(contacts, join_by(part_id)) %>% 
  mutate(count = replace_na(count, 0),
         date = as.Date(sday_id, format = "%Y.%m.%d"),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date)) %>% 
  group_by(fortnight) %>%
  summarise(mean_contacts = mean(count, na.rm = TRUE)) %>% 
  left_join(dates, by = "fortnight")

ggplot() +
  geom_line(data = combined, aes(x = date, y = mean_contacts), colour = "blue", lty = 2) +
  scale_x_date(date_breaks = "1 month") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) + 
  labs(x = "Date", y = "Mean Number of Contacts")

## age-stratified example
combined <- survey %>% 
  left_join(participants, join_by(part_id)) %>% 
  left_join(contacts, join_by(part_id)) %>% 
  mutate(count = replace_na(count, 0),
         date = as.Date(sday_id, format = "%Y.%m.%d"),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date)) %>% 
  mutate(agegp = case_when(part_age %in% c("Under 1", "0-4") ~ "0 to 4",
                           part_age %in% c("5-11", "12-17") ~ "5 to 14", #"12-17"
                           part_age %in% c("12-17", "18-19", "18-29", "25-34", "30-39", "35-44", "40-49") ~ "15 to 44", #"40-49"
                           part_age %in% c("40-49", "45-54", "50-59", "60-69") ~ "45 to 64", #"60-69"
                           part_age %in% c("60-69", "70-120") ~ "65+")) %>% 
  group_by(fortnight, agegp) %>%
  summarise(mean_contacts = mean(count, na.rm = TRUE)) %>% 
  mutate(agegp = factor(agegp,
                           levels = c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+"))) %>% 
  drop_na(agegp) %>%
  left_join(dates, by = "fortnight")

ggplot() +
  geom_line(data = combined, aes(x = date, y = mean_contacts, colour = agegp)) +
  scale_colour_viridis_d(option = "H", end = 0.8) + 
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


# imputation/processing for missing participant data


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


# Raw CoMix data ----
library(qs)

contacts <- qread(here("inst", "data", "dt_comix_share", "contacts.qs"))
households <- qread(here("inst", "data", "dt_comix_share", "households.qs"))
part_min <- qread(here("inst", "data", "dt_comix_share", "part_min.qs"))
part <- qread(here("inst", "data", "dt_comix_share", "part.qs"))

dates <- data.frame(date = seq(from = as.Date("23-03-2020", format = "%d-%m-%Y"), to = as.Date("02-03-2022", format = "%d-%m-%Y"), by = "day")) %>% 
  mutate(fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date))

## Scotland-only CoMix

## isolate Scotland
## area_3_name is the region variable the CoMix cleaning pipeline uses (amy's
## github); part_id is reused across countries and can change region between
## waves, so link the two tables on part_wave_uid instead
part_scotland <- part[area_3_name == "Scotland"]
contacts_scotland <- contacts[part_wave_uid %in% part_scotland$part_wave_uid]

## trim contacts: the top 1% of participants report 37% of all contact rows
## (median 2, max 2800), so cap each participant-wave, keeping individually
## reported contacts and then home/work/school first
trim_n <- 50

contacts_scotland <- contacts_scotland[
  order(cnt_mass == "individual", cnt_home, cnt_work, cnt_school, decreasing = TRUE),
  nth_cnt := seq_len(.N), by = part_wave_uid][nth_cnt <= trim_n]

age_limits <- c(0, 5, 11, 15, 25, 35, 45, 55, 65)
age_labels <- c("[0,5)", "[5,11)", "[11,15)", "[15,25)", "[25,35)",
                "[35,45)", "[45,55)", "[55,65)", "[65,Inf)")

## socialmixr expects one row per participant, so part_wave_uid becomes part_id.
## Ages are handed over as part_age_exact plus estimated bounds rather than as a
## part_age string for clean() to re-parse. clean() regenerates (and clobbers)
## part_age_est_min/max whenever a part_age column is present, so part_age is
## dropped here and these columns are left to pass through untouched. Exact ages
## go in part_age_exact; participants who only gave a range keep that range in
## the bounds, so assign_age_groups() can impute or sample from it.
part_sm <- part_scotland %>%
  select(part_id = part_wave_uid,
         part_age, part_age_group, part_age_est_min, part_age_est_max,
         date, weekday) %>%
  mutate(part_age_exact = suppressWarnings(as.integer(part_age)),
         part_age_group = as.character(part_age_group),
         part_age_est_min = coalesce(as.numeric(part_age_est_min),
                                     as.numeric(sub("-.*$", "", part_age_group))),
         part_age_est_max = coalesce(as.numeric(part_age_est_max),
                                     as.numeric(sub("^.*-", "", part_age_group))),
         dayofweek = match(weekday, c("Sunday", "Monday", "Tuesday", "Wednesday",
                                      "Thursday", "Friday", "Saturday")) - 1L,
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2)))) %>%
  select(-part_age) %>%
  group_by(fortnight) %>%
  mutate(mid_date = min(date) + floor((max(date) - min(date))/2)) %>%
  ungroup()

contacts_sm <- contacts_scotland %>%
  select(part_id = part_wave_uid,
         cnt_age, cnt_age_group,
         cnt_age_est_min, cnt_age_est_max) %>%
  mutate(cnt_age = as.character(cnt_age),
         cnt_age_group = as.character(cnt_age_group),
         cnt_age_est_min = as.numeric(cnt_age_est_min),
         cnt_age_est_max = as.numeric(cnt_age_est_max),
         cnt_age_exact = NA_integer_)

## Scottish population for survey.pop (the socialmixr default weights to the
## whole UK, via the country column)
pop_raw <- read_excel(here("inst", "data", "mid-year-population-estimates-time-series-data.xlsx"),
                      sheet = "Table 1", skip = 5)

scot_pop <- pop_raw %>%
  filter(`Area name` == "Scotland", Sex == "Persons", Year %in% c(2020, 2021, 2022)) %>% 
  group_by(`Area name`, `Sex`) %>% 
  summarise(across(`0`:`90 and over`, mean)) %>% 
  ungroup() %>% 
  select(`0`:`90 and over`) %>%
  mutate(across(everything(), as.numeric)) %>%
  pivot_longer(everything(), names_to = "age", values_to = "population") %>%
  mutate(age = as.numeric(sub(" and over", "", age))) %>% 
  arrange(age)

survey_pop <- data.frame(age = limits_to_age_groups(scot_pop$age, notation = "brackets"),
                         population = scot_pop$population)

# Fortnightly mean number of contacts for Scotland CoMix (total and by age group) ----
n_contacts <- contacts_sm %>%
  group_by(part_id) %>%
  summarise(n_contacts = n())

combined <- part_sm %>%
  left_join(n_contacts, join_by(part_id)) %>%
  mutate(n_contacts = replace_na(n_contacts, 0)) %>% 
  filter(!is.na(part_age_group)) %>% 
  mutate(agegp = case_when(part_age_group %in% c("0-4") ~ "0 to 4",
                           part_age_group %in% c("5-11", "12-17") ~ "5 to 14", #"12-17"
                           part_age_group %in% c("12-17", "18-29", "30-39", "40-49") ~ "15 to 44", #"40-49"
                           part_age_group %in% c("40-49", "50-59", "60-69") ~ "45 to 64", #"60-69"
                           part_age_group %in% c("60-69", "70-120") ~ "65+")) %>% 
  mutate(agegp = factor(agegp,
                        levels = c("0 to 4", "5 to 14", "15 to 44", "45 to 64", "65+")))

mean_total_scotland <- combined %>%
  group_by(fortnight, mid_date) %>%
  summarise(mean_contacts = mean(n_contacts), n_part = n()) %>%
  ungroup() %>%
  arrange(mid_date) %>% 
  left_join(dates, by = "fortnight")

mean_age_scotland <- combined %>%
  group_by(fortnight, mid_date, agegp) %>%
  summarise(mean_contacts = mean(n_contacts), n_part = n()) %>%
  ungroup() %>% 
  left_join(dates, by = "fortnight")

ggplot() +
  geom_line(data = mean_total_scotland, aes(x = date, y = mean_contacts),
            colour = "blue", lty = 2) +
  scale_x_date(date_breaks = "2 months") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "Date", y = "Mean Number of Contacts")

ggplot() +
  geom_line(data = mean_age_scotland, aes(x = date, y = mean_contacts, colour = agegp)) +
  scale_colour_viridis_d(option = "H", end = 0.8) +
  scale_x_date(date_breaks = "1 month") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) + 
  labs(x = "Date", y = "Mean Number of Contacts", colour = "Age Group")

# Fortnightly contact matrices for Scotland CoMix ----
part_split <- split(part_sm, part_sm$mid_date)

fortnight_survey_scotland <- lapply(part_split, function(p) {
  as_contact_survey(list(
    participants = as.data.frame(p),
    contacts = contacts_sm %>% filter(part_id %in% p$part_id) %>% as.data.frame()
  ))
})

fortnight_matrix_scotland <- lapply(fortnight_survey_scotland, function(s) {
  m <- s %>%
    assign_age_groups(age_limits = age_limits,
                      estimated_participant_age = "mean",
                      estimated_contact_age = "mean",
                      missing_participant_age = "remove",
                      missing_contact_age = "remove") %>%
    weigh_by_dayofweek() %>%
    compute_matrix()
  ## a handful of early/sparse fortnights have zero child participants,
  ## leaving NA rows/cols; symmetrise() refuses those, so leave them
  ## unsymmetrised here -- they get POLYMOD-patched below (see
  ## na_fortnight_scotland) and can be symmetrised again after patching
  if (anyNA(m$matrix)) {
    m
  } else {
    symmetrise(m, survey_pop = align_ages(survey_pop, m))
  }
})

## identify which matrices contain missing data: the first three fortnights have
## no child participants at all, so they need the same POLYMOD scaling applied to
## the UK matrices above
na_fortnight_scotland <- vapply(fortnight_matrix_scotland,
                                function(x) anyNA(x$matrix), logical(1))
names(which(na_fortnight_scotland))

## replace matrices with NA with closest complete matrix
fortnight_matrix_scotland <- fortnight_matrix_scotland[-54] # not date of interest
fortnight_matrix_scotland <- fortnight_matrix_scotland[-53] # not date of interest
fortnight_matrix_scotland[[1]][["matrix"]] <- fortnight_matrix_scotland[[4]][["matrix"]]
fortnight_matrix_scotland[[2]][["matrix"]] <- fortnight_matrix_scotland[[4]][["matrix"]]
fortnight_matrix_scotland[[3]][["matrix"]] <- fortnight_matrix_scotland[[4]][["matrix"]]
fortnight_matrix_scotland[[9]][["matrix"]] <- fortnight_matrix_scotland[[10]][["matrix"]] # replace with 10 because in same phase (phase 3 starts July)
fortnight_matrix_scotland[[52]][["matrix"]] <- fortnight_matrix_scotland[[51]][["matrix"]]

## example plot of single fortnightly matrix
fortnight_matrix_scotland[["2021-03-21"]]$matrix %>%
  as.data.frame() %>%
  rownames_to_column("part_age_group") %>%
  pivot_longer(-part_age_group, names_to = "contact_age_group", values_to = "mean_contacts") %>%
  mutate(part_age_group = factor(part_age_group, levels = age_labels),
         contact_age_group = factor(contact_age_group, levels = age_labels)) %>%
  ggplot(aes(x = part_age_group, y = contact_age_group, fill = mean_contacts)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D") +
  theme_bw() +
  labs(x = "Participant Age Group", y = "Contact Age Group", fill = "Mean Contacts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# Compare mean number of contacts (UK vs Scotland-only) ----
mean_age_all <- mean_age_scotland %>% 
  select(-c(mid_date, n_part)) %>% 
  mutate(source = "Scotland") %>% 
  bind_rows(combined %>% 
              mutate(source = "UK"))

ggplot(data = mean_age_all) +
  geom_line(aes(x = date, y = mean_contacts, colour = agegp, linetype = source)) +
  scale_colour_viridis_d(option = "H", end = 0.8) + 
  scale_x_date(date_breaks = "2 month") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) + 
  labs(x = "Date", y = "Mean Number of Contacts", colour = "Age Group", linetype = "Source") +
  facet_wrap(~agegp)
