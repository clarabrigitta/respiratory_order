# Build the fixed model inputs once, to run locally or for HPC ----
#
# Can run and keep objects in global environment or write a single .rds.
# For HPC: copy inst/outdata/hpc_inputs.rds to the cluster.

source(here::here("R", "library_and_scripts.r"))

# contact matrices (from explore_contacts.r) ----

age_limits <- c(0, 5, 11, 15, 25, 35, 45, 55, 65)

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
                                          age.limits = age_limits)
}

## fortnights 1-3 and 52 have missing participant data: carry over neighbours
n_fn <- length(fortnight_matrix)
fortnight_matrix[[1]]$matrix <- fortnight_matrix[[4]]$matrix
fortnight_matrix[[2]]$matrix <- fortnight_matrix[[4]]$matrix
fortnight_matrix[[3]]$matrix <- fortnight_matrix[[4]]$matrix
fortnight_matrix[[n_fn]]$matrix <- fortnight_matrix[[n_fn - 1]]$matrix

stopifnot(!any(vapply(fortnight_matrix, function(m) anyNA(m$matrix), logical(1))))

# population and births (from explore_scotland.r) ----

scot_population <- read_excel(here("inst", "data", "mid-year-population-estimates-time-series-data.xlsx"),
                              sheet = "Table 1", skip = 5) %>%
  filter(`Area name` == "Scotland", Sex == "Persons", Year %in% c(2020:2022)) %>%
  select(-`All Ages`) %>%
  pivot_longer(cols = 5:95, names_to = "age", values_to = "n") %>%
  mutate(agegp = case_when(age %in% as.character(0:4)   ~ "0-4",
                           age %in% as.character(5:10)  ~ "5-10",
                           age %in% as.character(11:14) ~ "11-14",
                           age %in% as.character(15:24) ~ "15-24",
                           age %in% as.character(25:34) ~ "25-34",
                           age %in% as.character(35:44) ~ "35-44",
                           age %in% as.character(45:54) ~ "45-54",
                           age %in% as.character(55:64) ~ "55-64",
                           age %in% c(as.character(65:89), "90 and over") ~ "65+")) %>%
  group_by(agegp, Year) %>%
  summarise(sum_n = sum(n), .groups = "drop") %>%
  group_by(agegp) %>%
  summarise(average_n = floor(mean(sum_n))) %>%  # average across 2020-2022
  ungroup() %>%
  mutate(agegp = factor(agegp, levels = c("0-4", "5-10", "11-14", "15-24", "25-34",
                                          "35-44", "45-54", "55-64", "65+"))) %>%
  arrange(agegp)

scot_births <- read_excel(here("inst", "data", "weekly-births-26-week-2.xlsx"),
                          sheet = "Table_1", skip = 4) %>%
  mutate(`Week beginning` = as.Date(`Week beginning`)) %>%
  filter(`Registration year` %in% c(2020:2022)) %>%
  rename(year = `Registration year`) %>%
  group_by(year) %>%
  summarise(births_annual = sum(`Births registered`)) %>%
  mutate(births_daily = births_annual / 365)

# save for hpc ----

dir.create(here("inst", "outdata"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(fortnight_matrix = fortnight_matrix,
             scot_population  = scot_population,
             scot_births      = scot_births,
             age_limits       = age_limits),
        file = here("inst", "outdata", "hpc_inputs.rds"))

cat("wrote inst/outdata/hpc_inputs.rds:",
    length(fortnight_matrix), "fortnightly contact matrices,",
    nrow(scot_population), "age groups\n")
