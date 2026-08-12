# modelled period ----
date_start <- as.Date("2020-03-23") # same window as used in explore_scotland.r
date_end   <- as.Date("2022-03-02")

# number of swabs tested each week (denominator for positivity) ----
swabs <- read_csv("inst/data/rcgp/Number of swabs by week.csv") %>%
  mutate(WeekBeginning = as.Date(ISOWeekStartDate, format = "%d/%m/%Y")) %>%
  rename(NumberSwabsPerWeek = `Number of swabs`) %>%
  select(WeekBeginning, ISOWeek, ISOYear, NumberSwabsPerWeek)

# load RCGP sentinel surveillance data by virus ----
## PositivityThisWeek" in the "percentage of positive samples" files is actually a count of positives,
## not a percentage (confirmed against "Number of samples by virus by week.csv" and the metadata sheet)
data <- read_csv("inst/data/rcgp/Number of samples by virus by week.csv") %>%
  mutate(WeekBeginning = as.Date(`Week commencing`, format = "%d/%m/%Y")) %>%
  rename(NumberPositivesPerWeek = NumberOfPositivesThisWeek) %>%
  filter(!Virus %in% c("None detected (negatives)", "AwaitingResult")) %>% # drop non-pathogen result categories
  select(WeekBeginning, Virus, NumberPositivesPerWeek) %>%
  arrange(Virus, WeekBeginning) %>%
  complete(WeekBeginning = seq(min(WeekBeginning), max(WeekBeginning), by = "week"),
           Virus,
           fill = list(NumberPositivesPerWeek = 0)) %>% # weeks with no positives are absent rather than zero
  left_join(swabs %>% select(WeekBeginning, NumberSwabsPerWeek), by = "WeekBeginning") %>%
  mutate(Positivity = 100 * NumberPositivesPerWeek / NumberSwabsPerWeek) %>%
  filter(WeekBeginning >= date_start, WeekBeginning <= date_end)

# plot RCGP data by virus ----
## ggplot - counts
ggplot(data) +
  geom_line(aes(x = WeekBeginning, y = NumberPositivesPerWeek, colour = Virus)) +
  theme_bw() +
  scale_colour_viridis_d(option = "H",
                         begin = 0,
                         end = 1) +
  scale_x_date(date_breaks = "3 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Number of Positives per Week")

## ggplot - positivity
## note: swab numbers collapse during this period (36-621 per week, median 108), so positivity is noisy
ggplot(data) +
  geom_line(aes(x = WeekBeginning, y = Positivity, colour = Virus)) +
  theme_bw() +
  scale_colour_viridis_d(option = "H",
                         begin = 0,
                         end = 1) +
  scale_x_date(date_breaks = "3 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  labs(x = "Week Beginning", y = "Positivity (%)")

# load RCGP data by virus and age ----
## virology-age-band.xlsx supersedes the per-age csv files: it covers 2019-10-07 to 2026-07-20 (so it
## spans the modelled period) and, unlike the csvs, carries the "None detected" rows, which give an
## age-specific denominator
## one sheet per age band, plus an "All ages" sheet (the four bands sum exactly to it)
## sheets have two rows of applied-filter text above the header, hence skip = 2
age_sheets <- c("Under 5" = "lt 5",
                "5-18"    = "between 5 & 18",
                "19-64"   = "between 19 & 64",
                "65+"     = "gt 64")

data_age_raw <- bind_rows(lapply(names(age_sheets), function(a) {
  read_excel("inst/data/rcgp/virology-age-band.xlsx", sheet = age_sheets[[a]], skip = 2) %>%
    mutate(AgeGroup = a)})) %>%
  mutate(WeekBeginning = as.Date(`Week commencing`),
         Virus = sub("\\*+$", "", VirusDisplayName)) %>% # strip footnote markers, e.g. "Other Coronavirus**"
  rename(NumberPositivesPerWeek = NumberOfPositivesThisWeek) %>%
  select(WeekBeginning, AgeGroup, Virus, NumberPositivesPerWeek)

## age-specific denominator: all results for that age band and week, i.e. positives plus "None detected"
## note: this slightly exceeds the number of samples where one sample has multiple detections (against the
## all-age swab counts the equivalent total runs ~2% high, median ratio 1.02), so treat it as a close proxy
samples_age <- data_age_raw %>%
  group_by(WeekBeginning, AgeGroup) %>%
  summarise(NumberSamplesPerWeek = sum(NumberPositivesPerWeek), .groups = "drop")

data_age <- data_age_raw %>%
  filter(Virus != "None detected") %>% # drop the negatives now they have been used as denominator
  arrange(AgeGroup, Virus, WeekBeginning) %>%
  complete(WeekBeginning = seq(min(WeekBeginning), max(WeekBeginning), by = "week"),
           AgeGroup,
           Virus,
           fill = list(NumberPositivesPerWeek = 0)) %>% # weeks with no positives are absent rather than zero
  left_join(samples_age, by = c("WeekBeginning", "AgeGroup")) %>%
  mutate(AgeGroup = factor(AgeGroup, levels = c("Under 5", "5-18", "19-64", "65+")), # set after the join, which would coerce a factor back to character
         Positivity = 100 * NumberPositivesPerWeek / NumberSamplesPerWeek) %>% # NA where no samples taken
  filter(WeekBeginning >= date_start, WeekBeginning <= date_end) # drop this line for the full 2019-2026 extent

# plot RCGP data by virus and age ----
## ggplot - counts
ggplot(data_age) +
  geom_line(aes(x = WeekBeginning, y = NumberPositivesPerWeek, colour = Virus)) +
  theme_bw() +
  scale_colour_viridis_d(option = "H",
                         begin = 0,
                         end = 1) +
  scale_x_date(date_breaks = "3 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  facet_wrap(~AgeGroup, ncol = 1, scales = "free_y") +
  labs(x = "Week Beginning", y = "Number of Positives per Week", colour = "Virus")

## ggplot - positivity
## note: age-band denominators are very small over this period (median samples per week: 12 under 5,
## 8 for 5-18, 57 for 19-64, 22 for 65+), so the paediatric bands in particular swing to 0/100% on
## single samples - read these alongside the counts above
ggplot(data_age) +
  geom_line(aes(x = WeekBeginning, y = Positivity, colour = Virus)) +
  theme_bw() +
  scale_colour_viridis_d(option = "H",
                         begin = 0,
                         end = 1) +
  scale_x_date(date_breaks = "3 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  facet_wrap(~AgeGroup, ncol = 1, scales = "free_y") +
  labs(x = "Week Beginning", y = "Positivity within Age Group (%)", colour = "Virus")
