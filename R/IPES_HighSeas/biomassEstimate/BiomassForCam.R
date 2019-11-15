#####################################################################################################
#
# BiomassForCam.R
# 2019-11-08
# Written by Jackie King and Jennifer Boldt
# Updated by Erika Anderson
#
# R Code from biomass code chunks in IPES report
#
#
#####################################################
# libraries for IPES_REPORT
# likely not all are necessary for this script

# load R libraries
library(csasdown) # to build document
library(plyr) # used in biomass estimate
library(tidyverse) # data manipulation and graphing
library(magrittr) # pipes
library(car) # statistics
library(RODBC) # database connection
library(kableExtra) # nice tables
library(lubridate) # functions with dates
library(readxl) # read excel files
library(stringr) # string manipulation
library(naniar) # function replace_with_na
library(cowplot) # multiple plots in one
library(float) # fix the table-page breaks
library(broom) # display linear model on LW graphs easily

################################ need up to date csv files to run ##################################

# based on original R code provided by JK from 2018 report
# original R code in file Input/2018/estimatingbiomassv5.R
# used inputs as csv files from MS Access queries
# net mensuration for 2019 not in main database, just EA's version

###########################

#from JK_VIEW_IPES_CPUE_BiomassEst
#only has tows where a catch occured, does not fully contain zero catches by species
fish <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                            "Biomass_Estimates.csv"))

# create simplier species code column
# create year and tow column for join to compare multiple years
fish %<>%
  dplyr::rename(SPECIES_CODE_orig = SPECIES_CODE) %>%
  mutate(SPECIES_CODE = as.integer(str_sub(SPECIES_CODE_orig, 1, 3))) %>%
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER))

##need to pull in all tows conducted from JK_VIEW_IPES_TRAWL_TOWS
#identify unique tows with year tow column
haul <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                             "Trawl_Tows.csv"))
haul %<>% 
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER))

#JB added
##JK: this is required because the actual volume sampled is not the volume of the block, but the volume of the tow!!!!  Derp.
volswept <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                                "VolumeSwept_ByStratum_DayNight.csv"), 
                     stringsAsFactors = FALSE)
haul <- merge(haul,volswept,by = c("TRIP_YEAR","STRATUM","DayNight"), all.x = TRUE)

##defining q
q_salmon <- 0.40 #assumed for salmon as per Volvenko. 2003. NPAFC Doc. 729. 32p.
q_herring <- 1.0 #otherwise q=1 for herring

##first summarize number of hauls per day/night/ headrope, stratum by year
summary_hauls <- haul %>%
  group_by(TRIP_YEAR, STRATUM, TARGET_HEADLINE_HEIGHT, DayNight) %>%
  summarise(strata_num = n()) 
annual_hauls <- haul %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_N = n())

##since in 2018 given poor weather there are several headrope, day, stratum combos that have only n=1 tow
##calculating biomass based on day and stratum only --- dropping headrope depth
##if you don't, then variance cannot be calculated for n = 1 tow



############## 
#EA: create biomass function to simplfy code 
#############
biomassFunction <- function(thisSpeciesCode, fish, haul, q_value) {
  
  # subset
  df <- fish %>% filter(SPECIES_CODE == thisSpeciesCode) %>%
    # create new column to use so function works on herring and salmon
    mutate(CatchWtCPUE_kg_km3 = 
             if_else(SPECIES_CODE == 96, 
                     TotalCatchWtCPUE_kg_km3, JuvCatchWtCPUE_kg_km3))
  
  
  #------------joining with the total number of tows, then adding in zeros where there is no cpue, filling in species code
  cpue_df <- left_join(haul, df, by = "yearTow", all.x = T, suffix = c("", ".fish"))
  cpue_df$CatchWtCPUE_kg_km3[is.na(cpue_df$CatchWtCPUE_kg_km3)] <- 0
  cpue_df$TotalCatchWtCPUE_kg_km3[is.na(cpue_df$TotalCatchWtCPUE_kg_km3)] <- 0
  cpue_df$SPECIES_CODE[is.na(cpue_df$SPECIES_CODE)] <- thisSpeciesCode
  cpue_df$CatchWtCPUE_kg_km3 <- cpue_df$CatchWtCPUE_kg_km3/1000  ##converting to tonnes
  
  ##JB added volume swept and put q here for mean, var, and sd
  ##JK: forgot to ask why, but q needs to be included in these three and not just applied at the end to the 
  ##JK: biomass estimation
  ##JK:  I originally had "mu_cpue_pink$strata_biomass = mu_cpue_pink$strata_cpue * mu_cpue_pink$strata_vol * q_salmon"
  
  #(A) calculate the biomass by year and day/night for each strata
  ## (A.1) summarize cpue and it's variance both accounting for catchability, count number of tows (n), and get the swept_volume (v) by strata
  mu_cpue_df <- cpue_df %>%
    group_by(TRIP_YEAR, STRATUM, DayNight) %>%
    summarise(strata_cpue = mean(CatchWtCPUE_kg_km3*q_value), 
              cpue_var = var(CatchWtCPUE_kg_km3*q_value),
              strata_num = n(),
              volswept = mean(SumOfOfficialVolumeSwept_km3)) %>%  ##taking the mean because in actual fact SumOfficialVolumeSwept is already the 
    ##the strata sum -- and since it's value is repeated for each tow, need to take
    ## the mean to get the strata sum back.
    ungroup()
  
  ## (A.2) assign the total volume (V) per stratum
  mu_cpue_df$strata_vol <- with(mu_cpue_df, ifelse(STRATUM == 504, 47.04, ##assigning stratum volume
                                                   ifelse(STRATUM == 505, 73.44, 
                                                          ifelse(STRATUM == 506,24.96,
                                                                 ifelse(STRATUM == 507, 40.32,
                                                                        ifelse(STRATUM == 508, 37.92,
                                                                               ifelse(STRATUM == 509, 44.64,
                                                                                      ifelse(STRATUM == 510, 91.20, 121.92))))))))
  ##(B) Calculate annual stratum biomass and variance
  ##(B.1) Biomass
  mu_cpue_df$strata_biomass = (mu_cpue_df$strata_cpue * mu_cpue_df$strata_vol)##this is where I had q previously 
  
  ##sample variance as per Thompson, S.K. 1992. Sampling. John Wiley and Sons, Inc. New York. 343 p.
  ##JK:  below was incorrect, because it assumed that the sampled volume was the block volume, when it is the tow volume
  #mu_cpue_pink$strata_N = mu_cpue_pink$strata_vol/(4*4*(30/1000))  #total number of 4x4 km blocks fished to 30 m per stratum 
  ##JB says:  use strata volume for Nh and volume swept for nh - except the last nh in the formula, and it should be number of tows.
  ##JK:  ie. Nh is strata volume as before with volume swept subtracted; then the variance is divided by the number of tows
  
  ##(B.2) Variance
  mu_cpue_df$strata_variance = mu_cpue_df$strata_vol * (mu_cpue_df$strata_vol - mu_cpue_df$volswept) *
    (mu_cpue_df$cpue_var / mu_cpue_df$strata_num)
  
  ## (C) sum up the biomass and variance estimates across all strata
  biomass_df <- mu_cpue_df %>%
    group_by(TRIP_YEAR, DayNight) %>%
    summarise(
      annual_biomass = sum(strata_biomass),
      annual_biomass_variance = sum(strata_variance),
      annual_biomass_num = sum(strata_num)
    ) %>%
    ungroup()
  
  ##(D) Then based on the annual biomass variance across all strata (in C above), calculate the standard deviation, SE, CV and CI            
  biomass_df$annual_biomass_sd <-
    sqrt(biomass_df$annual_biomass_variance)
  biomass_df$annual_biomass_cv <-
    biomass_df$annual_biomass_sd / biomass_df$annual_biomass
  biomass_df$annual_biomass_se <-
    biomass_df$annual_biomass_sd / (sqrt(biomass_df$annual_biomass_num))
  biomass_df$annual_biomass_LCI <-
    ifelse(((biomass_df$annual_biomass_sd * (-1.96) + biomass_df$annual_biomass)) < 0, 0,
           (biomass_df$annual_biomass_sd * (-1.96) + biomass_df$annual_biomass
           )) ##if LCI<0 then assign zero
  biomass_df$annual_biomass_UCI <-
    biomass_df$annual_biomass_sd * (1.96) + biomass_df$annual_biomass
  biomass_df$species_code <- thisSpeciesCode
  
  return(biomass_df)
  
}

#############

# apply biomass function to species
biomass_ck <- biomassFunction(124, fish, haul, q_salmon)
biomass_cm <- biomassFunction(112, fish, haul, q_salmon)
biomass_co <- biomassFunction(115, fish, haul, q_salmon)
biomass_pk <- biomassFunction(108, fish, haul, q_salmon)
biomass_se <- biomassFunction(118, fish, haul, q_salmon)
biomass_herring <- biomassFunction(96, fish, haul, q_herring)

# JB and JK only want daytime for salmon and nightime for herring
biomass_ck %<>% filter(DayNight == "Day")
biomass_cm %<>% filter(DayNight == "Day")
biomass_co %<>% filter(DayNight == "Day")
biomass_pk %<>% filter(DayNight == "Day")
biomass_se %<>% filter(DayNight == "Day")
biomass_herring %<>% filter(DayNight == "Night")

# bind species into one dataframe
biomass <- rbind(biomass_ck, biomass_cm, biomass_co, biomass_pk, biomass_se, biomass_herring)

# add common name and reorganize for table
biomass_display <- biomass %>%
  mutate(., species_code = if_else(species_code == 96, "096", as.character(species_code))) %>%
  left_join(., speciesCodes, by = c("species_code" = "SPECIES_CODE")) %>%
  mutate(COMMON_NAME = str_to_title(COMMON_NAME),
         annual_biomass = round(annual_biomass, 2),
         annual_biomass_cv = round(annual_biomass_cv, 2),
         annual_biomass_se = round(annual_biomass_se, 2),
         annual_biomass_LCI = round(annual_biomass_LCI, 2),
         annual_biomass_UCI = round(annual_biomass_UCI, 2)) %>%
  dplyr::rename(Species = COMMON_NAME,
                BiomassTonnes = annual_biomass,
                BiomassCV = annual_biomass_cv,
                BiomassSE = annual_biomass_se,
                LowerCI = annual_biomass_LCI ,
                UpperCI = annual_biomass_UCI) %>%
  filter(TRIP_YEAR == surveyYear) %>%
  select(Species, BiomassTonnes, BiomassCV, BiomassSE, LowerCI, UpperCI)

