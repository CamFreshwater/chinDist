### Simulate biomass estimates to quantify relative benefits of increasing 
# sampling effort within a strate

###### IMPORT DATA AND CLEAN
library(tidyverse)
library(ggplot2)
fish <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                            "Biomass_Estimates.csv"), stringsAsFactors = F) 

haul <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                            "Trawl_Tows.csv"), stringsAsFactors = F) %>%  
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER)) 

df <- fish %>% 
  dplyr::rename(SPECIES_CODE_orig = SPECIES_CODE) %>%
  mutate(SPECIES_CODE = as.integer(str_sub(SPECIES_CODE_orig, 1, 3))) %>%
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER),
         CatchWt = 
           if_else(SPECIES_CODE == 96, 
                   TotalCatchWtCPUE_kg_km3, JuvCatchWtCPUE_kg_km3) / 1000) %>% 
  select(TRIP_YEAR, STRATUM, yearTow, SPECIES_CODE, CatchWt)

# add zeros to catch data and merge with hauls
cpue_df <- haul %>% 
  select(yearTow, TRIP_YEAR, STRATUM, DayNight) %>% 
  left_join(., df %>% select(yearTow, SPECIES_CODE), by = "yearTow") %>% 
  expand(., nesting(yearTow, TRIP_YEAR, STRATUM, DayNight), SPECIES_CODE) %>% 
  filter(!is.na(SPECIES_CODE)) %>% 
  left_join(., df, by = c("yearTow", "TRIP_YEAR", "STRATUM", "SPECIES_CODE")) %>% 
  complete(yearTow, SPECIES_CODE, fill = list(CatchWt = 0)) %>% 
  mutate(nonZero = case_when(
    CatchWt > 0 ~ 1,
    CatchWt == 0 ~ 0
  )) %>% 
  group_by(SPECIES_CODE) %>% 
  nest() %>% 
  #don't account for year/strata now, but could as fixed or mixed effects
  mutate(m1 = map(data, ~glm(nonZero ~ 1, data = ., 
                             family = binomial(link = logit))),
         m2 = map(data, 
                  ~glm(CatchWt ~ 1, data = subset(., nonZero == 1), 
                       family = Gamma(link = log)))
         ) %>% 
  mutate(tidy1 = map(m1, broom::tidy),
         tidy2 = map(m2, broom::tidy))

cpue_df %>% 
  unnest(tidy1, tidy2) %>%
  select(species = SPECIES_CODE, data, binMu = estimate, binSig = std.error, 
         gamMu = estimate1, gamSig = std.error1) %>% 
  pivot_longer(-c(species, data), names_to = "parameter", 
               values_to = "estimate") %>% 
  mutate(model = case_when(
    grepl("gam", parameter) ~ "gamma",
    TRUE ~ "binomial"),
    parameter = case_when(
      grepl("Mu", parameter) ~ "mu",
      grepl("Sig", parameter) ~ "sigma"
    )) %>% 
  glimpse()

  grepl("PIT", stock) ~ "LWFR-Su"
  
relig_income
relig_income %>%
  pivot_longer(-religion,  names_to = "income", values_to = "count")


fish_encounters
fish_encounters %>%
  pivot_wider(names_from = station, values_from = seen)

#1. Use gamma hurdle models to estimate parameters describing data then simulate

#2. Sample from distribution at different levels
#3. Calculate CV following biomass extrapolotation


