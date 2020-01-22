## Multinomial model fit
# January 20, 2020
# Fit multinomial model to stock composition data 
# (assumes perfect GSI assignment)

library(tidyverse)
library(TMB)
library(ggplot2)


reg3_long <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "reg3RollUpCatchProb.RDS")) %>% 
  #remove all secondary IDs and uncertain IDs (less than 80%)
  group_by(flatFileID) %>% 
  mutate(maxProb = max(aggProb)) %>% 
  filter(!aggProb != maxProb, 
         !maxProb < 0.8) %>% 
  ungroup() %>% 
  # aggregate as necessary
  mutate(
    regName = case_when(
      regName %in% c("Columbia", "Snake") ~ "Columbia",
      regName %in% c("Coastal Washington", "Washington Coast", 
                     "Alaska South SE", "North/Central BC", "SOG") ~ "Other",
      TRUE ~ regName
    ),
    pres = 1) %>%
  dplyr::select(flatFileID, statArea, year, month, regName, pres)

gsi_wide <- reg3_long %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  filter(year == "2007")

obs_mat <- gsi_wide %>% 
  select()

compile("C:/github/juvenile-salmon-index/R/multinomialPractice/multinomial_generic.cpp")
dyn.load(dynlib("C:/github/juvenile-salmon-index/R/multinomialPractice/multinomial_generic"))




gsi_long_agg <- readRDS(here::here("data", "longGSI_reg4.rds")) %>% 
  filter(age == "J", 
         station_id %in% juv$station_id) %>% 
  mutate(jdayZ = as.vector(scale(jday)[,1]),
         present = 1,
         #consolidate northern aggregates because rel. rare
         agg = case_when(
           Region4Name %in% c("NBC", "SEAK", "CoastUS") ~ "Other",
           TRUE ~ Region4Name),
         season = as.factor(case_when(
           month %in% c("2", "3") ~ "winter",
           month %in% c("5", "6", "7", "8") ~ "summer",
           month %in% c("9", "10", "11" , "12") ~ "fall")),
         year = as.factor(year)
  ) %>% 
  select(-Region4Name, -ship_fl, -c(xUTM_start:age), -c(aggProb:maxProb)) %>%
  mutate(season = fct_relevel(season, "fall", after = 1)) %>% 
  left_join(., juv %>% select(station_id, week), by = "station_id") %>% 
  distinct()

year_aggs <- expand.grid(year = unique(gsi_long_agg$year), 
                         agg = unique(gsi_long_agg$agg), 
                         present = 1)

## Replace year/aggregate combinations with 0 catches with 1s and spread to wide
# format
gsi_wide <- full_join(gsi_long_agg, year_aggs, 
                      by = c("year", "agg", "present")) %>% 
  pivot_wider(., names_from = agg, values_from = present) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) 