### Simulate biomass estimates to quantify relative benefits of increasing 
# sampling effort within a strate

###### IMPORT DATA AND CLEAN
library(tidyverse)
library(ggplot2)
fish <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
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


df <- fish %>% 
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
  left_join(., df %>% select(yearTow, SPECIES_CODE), by = "yearTow") %>% 
  expand(., 
         nesting(yearTow, TRIP_YEAR, STRATUM, DayNight, meanSweepVol), 
         SPECIES_CODE) %>% 
  filter(!is.na(SPECIES_CODE)) %>% 
  left_join(., df, by = c("yearTow", "TRIP_YEAR", "STRATUM", "SPECIES_CODE")) %>% 
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
# %>% 
## REPLACED BY ABOVE
#   group_by(STRATUM) %>%
#   #take average of volume that was swept by tow ACROSS YEARS
#   mutate(volSweptMean = mean(volSweptStrata)) %>% 
#   select(-volSweptStrata) %>% 
#   ungroup() 
  

#1. Use gamma hurdle models to estimate parameters describing data then simulate
modfits <- cpue_df %>% 
  group_by(SPECIES_CODE) %>% 
  nest()  %>% 
  #don't account for year/strata now, but could as fixed or mixed effects
  mutate(m1 = map(data, ~glm(nonZero ~ STRATUM, data = ., 
                             family = binomial(link = logit))),
         m2 = map(data, 
                  ~glm(CatchWt ~ STRATUM, data = subset(., nonZero == 1), 
                       family = Gamma(link = log)))
         ) %>% 
  mutate(tidy1 = map(m1, broom::tidy),
         tidy2 = map(m2, broom::tidy),
         summ = map(m2, summary),
         gamDisp = map_dbl(summ, "dispersion"))


  
  
# Best practice would be to simulate data by specifying trials; each trial 
# generate an estimate of phi and beta (rate) by drawing from normal 
# distribution describing mu and se of parameter estimates. Then use these values
# to generate 1000 fishing events (i.e. 0/1 and estimate of catch if non-0). 
# Eventually this should account for strata and year effects as well.

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
M <- 20 #number of trials
N <- 200 #number of draws per trial (i.e. events per strata per trial)
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
  left_join(., volPars, by = "strata")

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
#NOTE: grouping of species, trials, years etc. is done external to function
calcBiomass <- function(sampledWt) {
  tt <- trialsOut %>% 
    filter(sampleSet == "20") %>%
    group_by(trial, strata) %>% #remove trial
    mutate(meanCPUE = mean(catchWt * q_value),
           varCPUE = var(catchWt * q_value)) %>% 
    ungroup() %>% 
    group_by(trial) %>% #remove
    mutate(strataBiomass = meanCPUE * strataVol,
           strataVar = strataVol * (strataVol - meanSweepVol) * 
             (varCPUE / sampleSet),
           annualBiomass = sum(strataBiomass),
           annualVariance = sum(strataVar),
           annualSamples = sum(sampleSet),
           annualBiomassSD = sqrt(annualVariance),
           annualBiomassCV = annualBiomassSD / annualBiomass,
           annualBiomassSE = annualBiomassSD / sqrt(annualSamples)) %>% 
    ungroup()
           
} 

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


ggplot(mu_cpue_df) + 
  geom_histogram(aes(x = volswept)) +
  facet_wrap(~STRATUM)

ggplot(tt) + 
  geom_histogram(aes(x = strataVar)) +
  facet_wrap(~strata)

mu_cpue_df %>% 
  select(volswept, STRATUM, TRIP_YEAR, DayNight) %>% 
  distinct() %>% 
  tail
  
#2. Sample from distribution at different levels
#3. Calculate CV following biomass extrapolotation


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