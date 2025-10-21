comix <- read_csv("data/contact_matrices_9_periods.csv") %>% 
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

comix_kids <- read_csv("data/contact_matrices_9_periods.csv") %>% 
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
  geom_line(data = contacts_daily, aes(x = date, y = mean_contacts)) +
  geom_line(data = data %>% filter(!Pathogen %in% c("Rhinovirus", "Adenovirus")), aes(x = WeekBeginning, y = NumberCasesPerWeek, colour = Pathogen)) +
  theme_bw()

ggplot() +
  geom_line(data = contacts_daily_kids, aes(x = date, y = contacts, colour = `Participant age`)) +
  scale_color_viridis(discrete = T, option = "D") +
  theme_bw()
