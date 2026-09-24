# Build the fixed model inputs once, to run locally or for HPC ----
#
# Can run and keep objects in global environment or write a single .rds.
# For HPC: copy inst/outdata/hpc_inputs.rds to the cluster.

source(here::here("R", "library_and_scripts.r"))

# contact matrices: Scotland-only CoMix subset (from explore_contacts.r) ----
#
# The UK-wide pipeline this replaces is kept, commented out, at the end of this
# section. Object names are unchanged (fortnight_survey, fortnight_matrix,
# age_limits), so nothing downstream needs editing -- only the contents differ.

age_limits <- c(0, 5, 11, 15, 25, 35, 45, 55, 65)

## the model window, expressed as the fortnight labels model_rcpp.r builds in
## `dates`. Scotland CoMix runs two waves past the window (Nov 2022), so the
## matrix list is keyed on these labels and trimmed to them, keeping the
## positional lookup in model_rcpp.r (dates$fortnight_n) aligned with the list
fit_start   <- as.Date("2020-03-23")
fit_end     <- as.Date("2022-03-02")
window_days <- seq(fit_start, fit_end, by = 1)
window_fortnights <- unique(paste(isoyear(window_days), "/",
                                  sprintf("%02d", ceiling(isoweek(window_days) / 2))))

## raw CoMix tables, read locally (replaces the zenodo get_survey() call)
part     <- qread(here("inst", "data", "dt_comix_share", "part.qs"))
contacts <- qread(here("inst", "data", "dt_comix_share", "contacts.qs"))

## isolate Scotland
## area_3_name is the region variable the CoMix cleaning pipeline uses (amy's
## github); part_id is reused across countries and can change region between
## waves, so link the two tables on part_wave_uid instead
part_scotland     <- part[area_3_name == "Scotland"]
contacts_scotland <- contacts[part_wave_uid %in% part_scotland$part_wave_uid]

## trim contacts: the top 1% of participants report 37% of all contact rows
## (median 2, max 2800), so cap each participant-wave, keeping individually
## reported contacts and then home/work/school first
trim_n <- 50

contacts_scotland <- contacts_scotland[
  order(cnt_mass == "individual", cnt_home, cnt_work, cnt_school, decreasing = TRUE),
  nth_cnt := seq_len(.N), by = part_wave_uid][nth_cnt <= trim_n]

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
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date) / 2)))) %>%
  select(-part_age)

contacts_sm <- contacts_scotland %>%
  select(part_id = part_wave_uid,
         cnt_age, cnt_age_group, cnt_age_est_min, cnt_age_est_max) %>%
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
  group_by(`Area name`, Sex) %>%
  summarise(across(`0`:`90 and over`, mean), .groups = "drop") %>%
  select(`0`:`90 and over`) %>%
  mutate(across(everything(), as.numeric)) %>%
  pivot_longer(everything(), names_to = "age", values_to = "population") %>%
  mutate(age = as.numeric(sub(" and over", "", age))) %>%
  arrange(age)

survey_pop <- data.frame(age = limits_to_age_groups(scot_pop$age, notation = "brackets"),
                         population = scot_pop$population)

## create surveys split by fortnight
## splitting on the fortnight label (rather than mid_date) names the list with
## the same labels dates$fortnight carries in model_rcpp.r, so the order is
## chronological and can be checked against the window below
part_split <- split(part_sm, part_sm$fortnight)

fortnight_survey <- lapply(part_split, function(p) {
  as_contact_survey(list(
    participants = as.data.frame(p),
    contacts = contacts_sm %>% filter(part_id %in% p$part_id) %>% as.data.frame()
  ))
})

## create contact matrix for each fortnight survey
fortnight_matrix <- lapply(fortnight_survey, function(s) {
  m <- s %>%
    assign_age_groups(age_limits = age_limits,
                      estimated_participant_age = "mean",
                      estimated_contact_age = "mean",
                      missing_participant_age = "remove",
                      missing_contact_age = "remove") %>%
    weigh_by_dayofweek() %>%
    compute_matrix()
  ## a handful of early/sparse fortnights have zero child participants, leaving
  ## NA rows/cols; symmetrise() refuses those, so leave them unsymmetrised here
  ## and patch them below
  if (anyNA(m$matrix)) m else symmetrise(m, survey_pop = align_ages(survey_pop, m))
})

## fortnights with no child participants at all: carry over the nearest complete
## neighbour. 2020/15 takes 2020/16 rather than 2020/14 because the two sit in
## the same phase (phase 3 starts July)
na_fortnight <- names(which(vapply(fortnight_matrix, function(m) anyNA(m$matrix), logical(1))))
cat("fortnights with missing contact data:", paste(na_fortnight, collapse = ", "), "\n")

patch <- c("2020 / 07" = "2020 / 10",
           "2020 / 08" = "2020 / 10",
           "2020 / 09" = "2020 / 10",
           "2020 / 15" = "2020 / 16",
           "2022 / 05" = "2022 / 04")

for (fn in names(patch)) {
  fortnight_matrix[[fn]]$matrix <- fortnight_matrix[[patch[[fn]]]]$matrix
}

## drop the fortnights outside the model window (the Nov 2022 CoMix waves)
stopifnot(all(window_fortnights %in% names(fortnight_matrix)))
fortnight_matrix <- fortnight_matrix[window_fortnights]

## model_rcpp.r indexes this list positionally through dates$fortnight_n, so it
## has to be exactly the window's fortnights, in chronological order
stopifnot(identical(names(fortnight_matrix), window_fortnights),
          !any(vapply(fortnight_matrix, function(m) anyNA(m$matrix), logical(1))))

## superseded: UK-wide CoMix pipeline ----
## uncomment this block (and comment out the Scotland block above) to rebuild
## hpc_inputs.rds from the UK-wide matrices instead
##
# ## load comix data from zenodo
# comix_uk <- get_survey("https://doi.org/10.5281/zenodo.4905745")
#
# ## load participant file to fill missing participant age
# participants <- read_csv("inst/data/CoMix_uk_participant_common.csv") %>% 
#   select(part_id, part_age) %>% 
#   mutate(part_age = recode(part_age, "Under 1" = "0-1")) %>% 
#   separate(part_age, into = c("part_age_est_min", "part_age_est_max"), sep = "-", convert = TRUE)
#
# ## fill in missing participant age
# comix_uk$participants <- comix_uk$participants %>%
#   left_join(participants, by = "part_id")
#
# ## define fortnight variable in participant data
# participants <- comix_uk$participants %>%
#   mutate(date = as.Date(sday_id, format = "%Y.%m.%d"),
#          fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2)))) %>% 
#   group_by(fortnight) %>% 
#   mutate(start_date = min(date),
#          end_date = max (date),
#          mid_date = start_date + floor((end_date - start_date)/2)) %>% 
#   ungroup() %>% 
#   group_split(mid_date) # separate data by fortnights
#
# ## define fortnight variable in contact data
# survey <- read_csv("inst/data/CoMix_uk_sday.csv") %>% 
#   select(part_id, sday_id)
#
# contacts <- comix_uk$contacts %>%
#   left_join(survey, join_by(part_id)) %>% 
#   mutate(date = as.Date(sday_id, format = "%Y.%m.%d"),
#          fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2)))) %>% 
#   group_by(fortnight) %>% 
#   mutate(start_date = min(date),
#          end_date = max (date),
#          mid_date = start_date + floor((end_date - start_date)/2)) %>% 
#   ungroup() %>% 
#   group_split(mid_date) # separate data by fortnights
#
# ## create surveys split by fortnight
# fortnight_survey <- list()
#
# for (i in 1:length(participants)) {
#   fortnight_survey[[i]] <- survey(participants[[i]], contacts[[i]])
# }
#
# ## create contact matrix for each fortnight survey
# fortnight_matrix <- list()
#
# for (i in 1:length(participants)) {
#   fortnight_matrix[[i]] <- contact_matrix(fortnight_survey[[i]],
#                                           age.limits = age_limits)
# }
#
# ## fortnights 1-3 and 52 have missing participant data: carry over neighbours
# n_fn <- length(fortnight_matrix)
# fortnight_matrix[[1]]$matrix <- fortnight_matrix[[4]]$matrix
# fortnight_matrix[[2]]$matrix <- fortnight_matrix[[4]]$matrix
# fortnight_matrix[[3]]$matrix <- fortnight_matrix[[4]]$matrix
# fortnight_matrix[[n_fn]]$matrix <- fortnight_matrix[[n_fn - 1]]$matrix
#
# stopifnot(!any(vapply(fortnight_matrix, function(m) anyNA(m$matrix), logical(1))))

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
