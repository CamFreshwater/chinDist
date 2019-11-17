### Simulate biomass estimates to quantify relative benefits of increasing 
# sampling effort within a strate

###### IMPORT DATA AND CLEAN
library(tidyverse)
library(ggplot2)
fishDat <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                            "Biomass_Estimates.csv"), stringsAsFactors = F) 

volswept <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                                "VolumeSwept_ByStratum_DayNight.csv"), 
                     stringsAsFactors = FALSE)
haulDat <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                            "Trawl_Tows.csv"), stringsAsFactors = F) %>%  
  merge(., volswept,by = c("TRIP_YEAR","STRATUM","DayNight"), 
        all.x = TRUE) %>% 
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER)) %>%
  #convoluted effort to calcuate average volume swept by strata since necessary 
  #to scale biomass estimates w/ sample size
  group_by(STRATUM, TRIP_YEAR, DayNight) %>%
  mutate(nTows = length(unique(TOW_NUMBER)),
         sweepVol = SumOfOfficialVolumeSwept_km3 / nTows) %>% 
  ungroup() %>% 
  group_by(STRATUM) %>% 
  mutate(meanSweepVol = mean(sweepVol)) %>% 
  ungroup()


fishDat <- fishDat %>% 
  dplyr::rename(SPECIES_CODE_orig = SPECIES_CODE) %>%
  mutate(SPECIES_CODE = as.integer(str_sub(SPECIES_CODE_orig, 1, 3))) %>%
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER),
         CatchWt = 
           if_else(SPECIES_CODE == 96, 
                   TotalCatchWtCPUE_kg_km3, JuvCatchWtCPUE_kg_km3) / 1000) %>% 
  select(TRIP_YEAR, STRATUM, yearTow, SPECIES_CODE, CatchWt)

# add zeros to catch data and merge with hauls
cpue_df <- haulDat %>% 
  select(yearTow, TRIP_YEAR, STRATUM, DayNight, meanSweepVol) %>% 
  left_join(., fishDat %>% select(yearTow, SPECIES_CODE), by = "yearTow") %>% 
  expand(., 
         nesting(yearTow, TRIP_YEAR, STRATUM, DayNight, meanSweepVol), 
         SPECIES_CODE) %>% 
  filter(!is.na(SPECIES_CODE)) %>% 
  left_join(., fishDat, by = c("yearTow", "TRIP_YEAR", "STRATUM", "SPECIES_CODE")) %>% 
  complete(yearTow, SPECIES_CODE, fill = list(CatchWt = 0)) %>% 
  #specify zero catches
  mutate(nonZero = 
           case_when(
             CatchWt > 0 ~ 1,
             CatchWt == 0 ~ 0
             ),
         #add catch ability coefficients
         q_value = 
           case_when(
             SPECIES_CODE == 96 ~ 1,
             TRUE ~ 0.4
           ),
         #add volume of strata
         strataVol = case_when(
           STRATUM == 504 ~ 47.04, 
           STRATUM == 505 ~ 73.44,
           STRATUM == 506 ~ 24.96,
           STRATUM == 507 ~ 40.32,
           STRATUM == 508 ~ 37.92,
           STRATUM == 509 ~ 44.64,
           STRATUM == 510 ~ 91.20, 
           STRATUM == 511 ~ 121.92)
         ) 

# Use gamma hurdle models to estimate parameters describing data then 
# simulate sampling events
simBiomass <- cpue_df %>% 
  split(., .$SPECIES_CODE) %>% 
  lapply(., function(x) {
    simTrials(x, M = 100, N = 1000, nVec = c(5, 10, 20, 50, 100))
    }) %>% 
  do.call(rbind.data.frame, .) %>%
  ungroup() %>% 
  mutate(nSetsPerStrata = as.factor(sampleSet),
         species = as.factor(species))

pdf(here::here("figs", "biomassSim", "annualCVEst.pdf"), height = 6, width = 8)
ggplot(simBiomass) +
  geom_boxplot(aes(x = nSetsPerStrata, y = annualBiomassCV)) +
  facet_wrap(~species, scales = "free_y")
dev.off()

# Function that fits the models and uses secondary helper functions to generate
# random data that is sampled
simTrials <- function(datIn, M = 100, N = 1000, 
                      nVec = c(5, 10, 20, 50, 100)) {
  #Fit models
  m1 <- glm(nonZero ~ STRATUM, data = datIn, family = binomial(link = logit))
  m2 <- glm(CatchWt ~ STRATUM, data = subset(datIn, nonZero == 1), 
            family = Gamma(link = log))
  
  strataLevels <- unique(datIn$STRATUM)
  newDat <- data.frame(STRATUM = strataLevels) #df of predictions
  phiStrata <- predict(m1, newdata = newDat, type='response') # on scale of response
  #gamma predictions
  muStrata <- predict(m2, newdata = newDat, type = 'response')
  #gamma distribution parameters
  gammaCoefs <- data.frame(strata = strataLevels,
                           shape = 1 / summary(m2)$dispersion) %>%
    mutate(rate = shape / muStrata)
 
  trialsOut  <- data.frame(strata = rep(strataLevels, each = M),
                           trial = rep(seq(1, M, by = 1), 
                                       times = length(strataLevels)),
                           #binomial coefficients first
                           phi = rep(phiStrata, each = M)) %>% 
    left_join(., gammaCoefs, by = "strata") %>% 
    #break into list because of issues expanding w/ random draws within a DF subset
    split(., .$trial) %>% 
    lapply(., function(x) {
      drawDist(x, N) %>% #generate samples based on model coefficients
        simSurvey(., nVec) #subsample based on input vector size 
    }) %>% 
    #recombine into DF
    do.call(rbind.data.frame, .) %>%
    ungroup() %>% 
    left_join(., volPars, by = "strata") %>% 
    select(-c(phi:rate)) 
  
  #calculate biomass based on CPUE
  biomassEst <- trialsOut %>% 
    calcBiomass() %>% 
    ungroup() %>% 
    mutate(species =  unique(datIn$SPECIES_CODE)) %>%  #add species ID
    select(species, sampleSet:annualBiomassSE) #reorder
  return(biomassEst)
}

# function that can be passed to lapply to generate draws from both distributions
# within a given trial
drawDist <- function(coefsDat, N) {
  draw <- seq(1, N, by = 1)
  
  # Sample binomial distribution within strata/trial N draw times
  tempBinomial <- coefsDat %>% 
    expand(nesting(strata, trial, phi, shape, rate), draw) %>% 
    group_by(strata) %>% 
    #don't account for year/strata now, but could as fixed or mixed effects
    mutate(nonZero = rbinom(N, size = 1, prob = phi)) 
  
  # Sample gamma distribution within strata/trial for nonzero trials within 
  # a strata
  tempGamma <- tempBinomial %>% 
    filter(nonZero == "1") %>% 
    group_by(strata) %>%
    mutate(towsWCatch = length(nonZero)) %>% 
    mutate(catchWt = rgamma(n =  towsWCatch, shape = shape, rate = rate)) %>% 
    select(strata, draw, catchWt)
  
  # Merge dataframes and export
  catchOut <- tempBinomial %>% 
    left_join(., tempGamma, by = c("strata", "draw")) %>% 
    replace_na(., list(catchWt = 0)) %>% 
    ungroup()
  return(catchOut)
}

#function that can be pass to lapply to simulate surveying across each trial 
#with different sample sizes 
simSurvey <- function(trialsDat, nVec) {
  simList <- vector(mode = "list", length = length(nVec))
  for (i in seq_along(simList)) {
    simList[[i]] <- trialsDat %>% 
      group_by(strata) %>% 
      sample_n(., size = nVec[i], replace = F) %>% 
      mutate(sampleSet = nVec[i]) %>% 
      ungroup()
  }
  do.call(rbind, simList)
}

#function modified from BiomassForCam.R to estimate mean and CV of biomass
calcBiomass <- function(sampledCatch) {
  sampledCatch %>% 
    # filter(sampleSet == "20", trial == "1") %>%
    group_by(sampleSet, trial, strata, strataVol, meanSweepVol) %>% 
    summarise(meanCPUE = mean(catchWt * q_value),
              varCPUE = var(catchWt * q_value)) %>% 
    mutate(strataBiomass = meanCPUE * strataVol,
           strataVar = strataVol * (strataVol - meanSweepVol) *
             (varCPUE / sampleSet)
           ) %>% 
    group_by(sampleSet, trial) %>%
    summarise(annualBiomass = sum(strataBiomass),
           annualVariance = sum(strataVar),
           annualSamples = sum(sampleSet),
           annualBiomassSD = sqrt(annualVariance),
           annualBiomassCV = annualBiomassSD / annualBiomass,
           annualBiomassSE = annualBiomassSD / sqrt(annualSamples)) 
} 

## pvalues for glms
pvals <- modfits %>% 
  mutate(mutate(tidy1 = map(m1, broom::tidy),
                tidy2 = map(m2, broom::tidy))) %>% 
  unnest(tidy1, tidy2) %>% 
  select(species = SPECIES_CODE, term, binMu = estimate,  
         binP = p.value, gamMu = estimate1, gamP = p.value1) %>% 
  pivot_longer(-c(species, term), names_to = "outputValue", 
               values_to = "estimate") %>%
  mutate(
    model = case_when(
      grepl("gam", outputValue) ~ "gamma",
      TRUE ~ "binomial"),
    value = case_when(
      grepl("Mu", outputValue) ~ "mu",
      grepl("P", outputValue) ~ "p-value"
    )) %>%
  select(-outputValue)

pvals %>% 
  filter(term == "STRATUM", 
         value == "p-value")
#generally strata effects are significant for the binomial model, but less 
#consistent for gamma

## Preliminary biomass sampling and calcs with one dataframe (i.e. species)

df <- data_frame(one = rep("hey", 10), two = seq(1:10), etc = "etc")

list_df <- list(df, df, df, df, df)
dfnames <- c("first", "second", "third", "fourth", "fifth")

dfs <- list_df %>% map2_df(dfnames,~mutate(.x, name=.y))

# simpler example with one species  
volPars <- cpue_df %>% 
  filter(SPECIES_CODE == "124") %>% 
  select(strata = STRATUM, q_value, strataVol, meanSweepVol) %>% 
  distinct()

dum <- cpue_df %>% 
  filter(SPECIES_CODE == "124")
m1 <- glm(nonZero ~ STRATUM, data = dum, family = binomial(link = logit))
m2 <- glm(CatchWt ~ STRATUM, data = subset(dum, nonZero == 1), 
          family = Gamma(link = log))

newDat <- data.frame(STRATUM = unique(dum$STRATUM))
#binomial predictions
phiStrata <- predict(m1, newdata=newDat, type='response') # on scale of response
#gamma predictions
muStrata <- predict(m2, newdata = newDat, type = 'response')
#gamma distribution parameters
gammaCoefs <- data.frame(strata = unique(dum$STRATUM),
                         shape = 1 / summary(m2)$dispersion) %>%
  mutate(rate = shape / muStrata)

# generate random draws for the binomial catch
M <- 100 #number of trials
N <- 1000 #number of draws per trial (i.e. events per strata per trial)
nVec <- c(5, 10, 20, 50, 100) #vector of survey sets (per strata)
trialsOut  <- data.frame(strata = rep(unique(dum$STRATUM), each = M),
                         trial = rep(seq(1, M, by = 1), 
                                     times = length(unique(dum$STRATUM))),
                         #binomial coefficients first
                         phi = rep(phiStrata, each = M)) %>% 
  left_join(., gammaCoefs, by = "strata") %>% 
  #break into list because of issues expanding w/ random draws within a DF subset
  split(., .$trial) %>% 
  lapply(., function(x) {
    drawDist(x, N) %>% #generate samples based on model coefficients
      simSurvey(., nVec) #subsample based on input vector size 
  }) %>% 
  #recombine into DF
  do.call(rbind.data.frame, .) %>%
  ungroup() %>% 
  left_join(., volPars, by = "strata") %>% 
  select(-c(phi:rate)) 

biomassEst <- trialsOut %>% 
  calcBiomass()

ggplot(biomassEst) + 
  geom_boxplot(aes(x = as.factor(sampleSet), y = annualVariance))