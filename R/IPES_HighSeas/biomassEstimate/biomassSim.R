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
                   TotalCatchWtCPUE_kg_km3, JuvCatchWtCPUE_kg_km3)) %>% 
  select(TRIP_YEAR, STRATUM, yearTow, SPECIES_CODE, CatchWt)
  
## Merge with haul
cpue_df <- left_join(haul, df, by = c("yearTow", "TRIP_YEAR", "STRATUM")) %>% 
  filter(SPECIES_CODE == "96") %>%
  complete(nesting(TOW_NUMBER, yearTow), SPECIES_CODE) %>% 
  # filter(!is.na(SPECIES_CODE)) %>%
  arrange(yearTow)
cpue_df$CatchWt[is.na(cpue_df$CatchWt)] <- 0
# cpue_df[is.na(cpue_df$SPECIES_CODE),] <- thisSpeciesCode
cpue_df$CatchWt <- cpue_df$CatchWt / 1000 

dum <- df %>% 
  filter(#STRATUM == "505",
         SPECIES_CODE == "96")
         #TRIP_YEAR == "2018")


#1. Identify plausible lognormal mean and sd of biomass for each species
ggplot(cpue_df) +
  geom_histogram(aes(x = CatchWt)) +
  facet_wrap(~SPECIES_CODE, scales = "free")

haul %>% 
  group_by(TRIP_YEAR, STRATUM) %>% 
  tally()

zeroinfl()

#2. Sample from distribution at different levels
#3. Calculate CV following biomass extrapolotation


