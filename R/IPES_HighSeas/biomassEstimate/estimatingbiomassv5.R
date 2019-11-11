## JB all species calculations of SDs and CIs (May 15, 2019)
## DR. ROOPER SAYS WE CAN'T SUM SDs OR SEs OR CIs  :( (May 15, 2019)
## King edited remaining, and cleaned up some code not used.


rm(list = ls(all=TRUE)); #Remove all the objects in the memory

setwd("C:/data/Documents/MS/Synoptic Pelagics Ecosystem Survey/2017&2018 Tech Report/Biomass Estimates")


######################
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)


#from JK_VIEW_IPES_CPUE_BiomassEst
#only has tows where a catch occured, does not fully contain zero catches by species
fish<- read.csv("Biomass_Estimates.csv", header=TRUE, sep=",")


##need to pull in all tows conducted from JK_VIEW_IPES_TRAWL_TOWS
#identify unique tows
haul <- read.csv("Trawl_Tows.csv", header=TRUE, sep = ",")
haul$haul.TOW_NUMBER <- haul$TOW_NUMBER

#JB added
##JK: this is required because the actual volume sampled is not the volume of the block, but the volume of the tow!!!!  Derp.
volswept<-read.csv("VolumeSwept_ByStratum_DayNight.csv", header=TRUE, sep=",")
haul<-merge(haul,volswept,by=c("TRIP_YEAR","STRATUM","DayNight"), all.x=TRUE)

#----------subsetting-----------------------------

pink <- fish %>% filter(SPECIES_CODE==108 ) 
chum <- fish %>% filter(SPECIES_CODE==112)  
coho <- fish %>% filter(SPECIES_CODE==115 ) 
sockeye <- fish %>% filter(SPECIES_CODE==118  )  
chinook <- fish %>% filter(SPECIES_CODE==124  ) 
herring <- fish %>% filter(SPECIES_CODE==96)    

#------------joining with the total number of tows, then adding in zeros where there is no cpue, filling in species code
cpue_pink <- left_join(haul, pink, by.x="haul.TOW_NUMBER",by.y="TOW_NUMBER",all.x=T)
cpue_pink$JuvCatchWtCPUE_kg_km3[is.na(cpue_pink$JuvCatchWtCPUE_kg_km3)]<-0
cpue_pink$TotalCatchWtCPUE_kg_km3[is.na(cpue_pink$TotalCatchWtCPUE_kg_km3)]<-0
cpue_pink$SPECIES_CODE[is.na(cpue_pink$SPECIES_CODE)]<-108
cpue_pink$JuvCatchWtCPUE_kg_km3 <- cpue_pink$JuvCatchWtCPUE_kg_km3/1000  ##converting to tonnes

cpue_chum <- left_join(haul, chum, by.x="haul.TOW_NUMBER",by.y="TOW_NUMBER",all.x=T)
cpue_chum$JuvCatchWtCPUE_kg_km3[is.na(cpue_chum$JuvCatchWtCPUE_kg_km3)]<-0
cpue_chum$TotalCatchWtCPUE_kg_km3[is.na(cpue_chum$TotalCatchWtCPUE_kg_km3)]<-0
cpue_chum$SPECIES_CODE[is.na(cpue_chum$SPECIES_CODE)]<-112
cpue_chum$JuvCatchWtCPUE_kg_km3 <- cpue_chum$JuvCatchWtCPUE_kg_km3/1000  ##converting to tonnes

cpue_coho <- left_join(haul, coho, by.x="haul.TOW_NUMBER",by.y="TOW_NUMBER",all.x=T)
cpue_coho$JuvCatchWtCPUE_kg_km3[is.na(cpue_coho$JuvCatchWtCPUE_kg_km3)]<-0
cpue_coho$TotalCatchWtCPUE_kg_km3[is.na(cpue_coho$TotalCatchWtCPUE_kg_km3)]<-0
cpue_coho$SPECIES_CODE[is.na(cpue_coho$SPECIES_CODE)]<-115
cpue_coho$JuvCatchWtCPUE_kg_km3 <- cpue_coho$JuvCatchWtCPUE_kg_km3/1000  ##converting to tonnes

cpue_sockeye <- left_join(haul, sockeye, by.x="haul.TOW_NUMBER",by.y="TOW_NUMBER",all.x=T)
cpue_sockeye$JuvCatchWtCPUE_kg_km3[is.na(cpue_sockeye$JuvCatchWtCPUE_kg_km3)]<-0
cpue_sockeye$TotalCatchWtCPUE_kg_km3[is.na(cpue_sockeye$TotalCatchWtCPUE_kg_km3)]<-0
cpue_sockeye$SPECIES_CODE[is.na(cpue_sockeye$SPECIES_CODE)]<-118
cpue_sockeye$JuvCatchWtCPUE_kg_km3 <- cpue_sockeye$JuvCatchWtCPUE_kg_km3/1000  ##converting to tonnes

cpue_chinook <- left_join(haul, chinook, by.x="haul.TOW_NUMBER",by.y="TOW_NUMBER",all.x=T)
cpue_chinook$JuvCatchWtCPUE_kg_km3[is.na(cpue_chinook$JuvCatchWtCPUE_kg_km3)]<-0
cpue_chinook$TotalCatchWtCPUE_kg_km3[is.na(cpue_chinook$TotalCatchWtCPUE_kg_km3)]<-0
cpue_chinook$SPECIES_CODE[is.na(cpue_chinook$SPECIES_CODE)]<-124
cpue_chinook$JuvCatchWtCPUE_kg_km3 <- cpue_chinook$JuvCatchWtCPUE_kg_km3/1000  ##converting to tonnes

cpue_herring <- left_join(haul, herring, by.x="haul.TOW_NUMBER",by.y="TOW_NUMBER",all.x=T)
cpue_herring$JuvCatchWtCPUE_kg_km3[is.na(cpue_herring$JuvCatchWtCPUE_kg_km3)]<-0
cpue_herring$TotalCatchWtCPUE_kg_km3[is.na(cpue_herring$TotalCatchWtCPUE_kg_km3)]<-0
cpue_herring$SPECIES_CODE[is.na(cpue_herring$SPECIES_CODE)]<-96
cpue_herring$TotalCatchWtCPUE_kg_km3 <- cpue_herring$TotalCatchWtCPUE_kg_km3/1000  ##converting to tonnes


##-----------------Estimating Biomass----------------------
##defining q
q_salmon <- 0.40 #assumed for salmon as per Volvenko. 2003. NPAFC Doc. 729. 32p.
q_herring <- 1.0 #otherwise q=1 for herring


##first summarize number of hauls per day/night/ headrope, stratum by year
summary_hauls <- haul %>%
  group_by(TRIP_YEAR, STRATUM, TARGET_HEADLINE_HEIGHT, DayNight) %>%
  summarise(strata_num=n()) 
annual_hauls <- haul %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_N=n())

##since in 2018 given poor weather there are several headrope, day, stratum combos that have only n=1 tow
##calculating biomass based on day and stratum only --- dropping headrope depth
##if you don't, then variance cannot be calculated for n=1 tow


##--------------pink
##JB added volume swept and put q here for mean, var, and sd
##JK: forgot to ask why, but q needs to be included in these three and not just applied at the end to the 
##JK: biomass estimation
##JK:  I originally had "mu_cpue_pink$strata_biomass = mu_cpue_pink$strata_cpue * mu_cpue_pink$strata_vol * q_salmon"

#(A) calculate the biomass by year and day/night for each strata
## (A.1) summarize cpue and it's variance both accounting for catchability, count number of tows (n), and get the swept_volume (v) by strata
mu_cpue_pink <- cpue_pink %>%
  group_by(TRIP_YEAR, STRATUM, DayNight) %>%
summarise(strata_cpue =mean(JuvCatchWtCPUE_kg_km3*q_salmon), 
          cpue_var =var(JuvCatchWtCPUE_kg_km3*q_salmon),
          strata_num=n(),
          volswept=mean(SumOfOfficialVolumeSwept_km3)) %>%  ##taking the mean because in actual fact SumOfficialVolumeSwept is already the 
  ##the strata sum -- and since it's value is repeated for each tow, need to take
  ## the mean to get the strata sum back.
  ungroup()
## (A.2) assign the total volume (V) per stratum
mu_cpue_pink$strata_vol <-with(mu_cpue_pink, ifelse(STRATUM==504, 47.04, ##assigning stratum volume
                                            ifelse(STRATUM==505, 73.44, 
                                            ifelse(STRATUM==506,24.96,
                                            ifelse(STRATUM==507, 40.32,
                                            ifelse(STRATUM==508, 37.92,
                                            ifelse(STRATUM==509, 44.64,
                                            ifelse(STRATUM==510, 91.20, 121.92))))))))
##(B) Calculate annual stratum biomass and variance
##(B.1) Biomass
mu_cpue_pink$strata_biomass = (mu_cpue_pink$strata_cpue * mu_cpue_pink$strata_vol)##this is where I had q previously 

##sample variance as per Thompson, S.K. 1992. Sampling. John Wiley and Sons, Inc. New York. 343 p.
##JK:  below was incorrect, because it assumed that the sampled volume was the block volume, when it is the tow volume
#mu_cpue_pink$strata_N = mu_cpue_pink$strata_vol/(4*4*(30/1000))  #total number of 4x4 km blocks fished to 30 m per stratum 
##JB says:  use strata volume for Nh and volume swept for nh - except the last nh in the formula, and it should be number of tows.
##JK:  ie. Nh is strata volume as before with volume swept subtracted; then the variance is divided by the number of tows

##(B.2) Variance
mu_cpue_pink$strata_variance = mu_cpue_pink$strata_vol*(mu_cpue_pink$strata_vol - mu_cpue_pink$volswept)*(mu_cpue_pink$cpue_var/mu_cpue_pink$strata_num)

## (C) sum up the biomass and variance estimates across all strata
biomass_pink <-mu_cpue_pink %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_biomass=sum(strata_biomass),
            annual_biomass_variance=sum(strata_variance),
            annual_biomass_num=sum(strata_num))

##(D) Then based on the annual biomass variance across all strata (in C above), calculate the standard deviation, SE, CV and CI            
biomass_pink$annual_biomass_sd<-sqrt(biomass_pink$annual_biomass_variance)
biomass_pink$annual_biomass_cv<-biomass_pink$annual_biomass_sd/biomass_pink$annual_biomass
biomass_pink$annual_biomass_se<-biomass_pink$annual_biomass_sd/(sqrt(biomass_pink$annual_biomass_num))
biomass_pink$annual_biomass_LCI<-ifelse(((biomass_pink$annual_biomass_sd*(-1.96)+biomass_pink$annual_biomass))<0,0,(biomass_pink$annual_biomass_sd*(-1.96)+biomass_pink$annual_biomass)) ##if LCI<0 then assign zero
biomass_pink$annual_biomass_UCI<-biomass_pink$annual_biomass_sd*(1.96)+biomass_pink$annual_biomass
biomass_pink$species_code <- 108

## (E) Ta da!  
write.table(biomass_pink, file="biomass_estimates_fromR.v3.csv", append=FALSE, row.names=FALSE, col.names=TRUE, sep=",") ##APPEND = FALSE TO CLEAR OUT AND START OVER

#--------------------------chum

#(A) calculate the biomass by year and day/night for each strata
## (A.1) summarize cpue and it's variance both accounting for catchability, count number of tows (n), and get the swept_volume (v) by strata
mu_cpue_chum <- cpue_chum %>%
  group_by(TRIP_YEAR, STRATUM, DayNight) %>%
  summarise(strata_cpue =mean(JuvCatchWtCPUE_kg_km3*q_salmon), 
            cpue_var =var(JuvCatchWtCPUE_kg_km3*q_salmon),
            strata_num=n(),
            volswept=mean(SumOfOfficialVolumeSwept_km3)) %>%  
  ungroup()

## (A.2) assign the total volume (V) per stratum
mu_cpue_chum$strata_vol <-with(mu_cpue_chum, ifelse(STRATUM==504, 47.04, 
                                                    ifelse(STRATUM==505, 73.44, 
                                                           ifelse(STRATUM==506,24.96,
                                                                  ifelse(STRATUM==507, 40.32,
                                                                         ifelse(STRATUM==508, 37.92,
                                                                                ifelse(STRATUM==509, 44.64,
                                                                                                                                                                              ifelse(STRATUM==510, 91.20, 121.92))))))))
##(B) Calculate annual stratum biomass and variance
##(B.1) Biomass
mu_cpue_chum$strata_biomass = (mu_cpue_chum$strata_cpue * mu_cpue_chum$strata_vol) 

## (B.2) Variance
mu_cpue_chum$strata_variance = mu_cpue_chum$strata_vol*(mu_cpue_chum$strata_vol - mu_cpue_chum$volswept)*(mu_cpue_chum$cpue_var/mu_cpue_chum$strata_num)

## (C) sum up the biomass and variance estimates across all strata
biomass_chum <-mu_cpue_chum %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_biomass=sum(strata_biomass),
            annual_biomass_variance=sum(strata_variance),
            annual_biomass_num=sum(strata_num))

##(D) Then based on the annual biomass variance across all strata (in C above), calculate the standard deviation, CI and CV
biomass_chum$annual_biomass_sd<-sqrt(biomass_chum$annual_biomass_variance)
biomass_chum$annual_biomass_cv<-biomass_chum$annual_biomass_sd/biomass_chum$annual_biomass
biomass_chum$annual_biomass_se<-biomass_chum$annual_biomass_sd/(sqrt(biomass_chum$annual_biomass_num))
biomass_chum$annual_biomass_LCI<-ifelse(((biomass_chum$annual_biomass_sd*(-1.96)+biomass_chum$annual_biomass))<0,0,(biomass_chum$annual_biomass_sd*(-1.96)+biomass_chum$annual_biomass)) ##if LCI<0 then assign zero
biomass_chum$annual_biomass_UCI<-biomass_chum$annual_biomass_sd*(1.96)+biomass_chum$annual_biomass
biomass_chum$species_code <- 112

## (E) Ta da!
write.table(biomass_chum, file="biomass_estimates_fromR.v3.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",") 

#--------------------------coho

#(A) calculate the biomass by year and day/night for each strata
## (A.1) summarize cpue and it's variance both accounting for catchability, count number of tows (n), and get the swept_volume (v) by strata
mu_cpue_coho <- cpue_coho %>%
  group_by(TRIP_YEAR, STRATUM, DayNight) %>%
  summarise(strata_cpue =mean(JuvCatchWtCPUE_kg_km3*q_salmon), 
            cpue_var =var(JuvCatchWtCPUE_kg_km3*q_salmon),
            strata_num=n(),
            volswept=mean(SumOfOfficialVolumeSwept_km3)) %>%  
  ungroup()

## (A.2) assign the total volume (V) per stratum
mu_cpue_coho$strata_vol <-with(mu_cpue_coho, ifelse(STRATUM==504, 47.04, 
                                                    ifelse(STRATUM==505, 73.44, 
                                                           ifelse(STRATUM==506,24.96,
                                                                  ifelse(STRATUM==507, 40.32,
                                                                         ifelse(STRATUM==508, 37.92,
                                                                                ifelse(STRATUM==509, 44.64,
                                                                                       ifelse(STRATUM==510, 91.20, 121.92))))))))
##(B) Calculate annual stratum biomass and variance
##(B.1) Biomass
mu_cpue_coho$strata_biomass = (mu_cpue_coho$strata_cpue * mu_cpue_coho$strata_vol) 

## (B.2) Variance
mu_cpue_coho$strata_variance = mu_cpue_coho$strata_vol*(mu_cpue_coho$strata_vol - mu_cpue_coho$volswept)*(mu_cpue_coho$cpue_var/mu_cpue_coho$strata_num)

## (C) sum up the biomass and variance estimates across all strata
biomass_coho <-mu_cpue_coho %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_biomass=sum(strata_biomass),
            annual_biomass_variance=sum(strata_variance),
            annual_biomass_num=sum(strata_num))

##(D) Then based on the annual biomass variance across all strata (in C above), calculate the standard deviation, CI and CV
biomass_coho$annual_biomass_sd<-sqrt(biomass_coho$annual_biomass_variance)
biomass_coho$annual_biomass_cv<-biomass_coho$annual_biomass_sd/biomass_coho$annual_biomass
biomass_coho$annual_biomass_se<-biomass_coho$annual_biomass_sd/(sqrt(biomass_coho$annual_biomass_num))
biomass_coho$annual_biomass_LCI<-ifelse(((biomass_coho$annual_biomass_sd*(-1.96)+biomass_coho$annual_biomass))<0,0,(biomass_coho$annual_biomass_sd*(-1.96)+biomass_coho$annual_biomass)) ##if LCI<0 then assign zero
biomass_coho$annual_biomass_UCI<-biomass_coho$annual_biomass_sd*(1.96)+biomass_coho$annual_biomass
biomass_coho$species_code <- 115

## (E) Ta da!
write.table(biomass_coho, file="biomass_estimates_fromR.v3.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",") 


#--------------------------sockeye

#(A) calculate the biomass by year and day/night for each strata
## (A.1) summarize cpue and it's variance both accounting for catchability, count number of tows (n), and get the swept_volume (v) by strata
mu_cpue_sockeye <- cpue_sockeye %>%
  group_by(TRIP_YEAR, STRATUM, DayNight) %>%
  summarise(strata_cpue =mean(JuvCatchWtCPUE_kg_km3*q_salmon), 
            cpue_var =var(JuvCatchWtCPUE_kg_km3*q_salmon),
            strata_num=n(),
            volswept=mean(SumOfOfficialVolumeSwept_km3)) %>%  
  ungroup()

## (A.2) assign the total volume (V) per stratum
mu_cpue_sockeye$strata_vol <-with(mu_cpue_sockeye, ifelse(STRATUM==504, 47.04, 
                                                    ifelse(STRATUM==505, 73.44, 
                                                           ifelse(STRATUM==506,24.96,
                                                                  ifelse(STRATUM==507, 40.32,
                                                                         ifelse(STRATUM==508, 37.92,
                                                                                ifelse(STRATUM==509, 44.64,
                                                                                       ifelse(STRATUM==510, 91.20, 121.92))))))))
##(B) Calculate annual stratum biomass and variance
##(B.1) Biomass
mu_cpue_sockeye$strata_biomass = (mu_cpue_sockeye$strata_cpue * mu_cpue_sockeye$strata_vol) 

## (B.2) Variance
mu_cpue_sockeye$strata_variance = mu_cpue_sockeye$strata_vol*(mu_cpue_sockeye$strata_vol - mu_cpue_sockeye$volswept)*(mu_cpue_sockeye$cpue_var/mu_cpue_sockeye$strata_num)

## (C) sum up the biomass and variance estimates across all strata
biomass_sockeye <-mu_cpue_sockeye %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_biomass=sum(strata_biomass),
            annual_biomass_variance=sum(strata_variance),
            annual_biomass_num=sum(strata_num))

##(D) Then based on the annual biomass variance across all strata (in C above), calculate the standard deviation, CI and CV
biomass_sockeye$annual_biomass_sd<-sqrt(biomass_sockeye$annual_biomass_variance)
biomass_sockeye$annual_biomass_cv<-biomass_sockeye$annual_biomass_sd/biomass_sockeye$annual_biomass
biomass_sockeye$annual_biomass_se<-biomass_sockeye$annual_biomass_sd/(sqrt(biomass_sockeye$annual_biomass_num))
biomass_sockeye$annual_biomass_LCI<-ifelse(((biomass_sockeye$annual_biomass_sd*(-1.96)+biomass_sockeye$annual_biomass))<0,0,(biomass_sockeye$annual_biomass_sd*(-1.96)+biomass_sockeye$annual_biomass)) ##if LCI<0 then assign zero
biomass_sockeye$annual_biomass_UCI<-biomass_sockeye$annual_biomass_sd*(1.96)+biomass_sockeye$annual_biomass
biomass_sockeye$species_code <- 118

## (E) Ta da!
write.table(biomass_sockeye, file="biomass_estimates_fromR.v3.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",") 



#--------------------------chinook

#(A) calculate the biomass by year and day/night for each strata
## (A.1) summarize cpue and it's variance both accounting for catchability, count number of tows (n), and get the swept_volume (v) by strata
mu_cpue_chinook <- cpue_chinook %>%
  group_by(TRIP_YEAR, STRATUM, DayNight) %>%
  summarise(strata_cpue =mean(JuvCatchWtCPUE_kg_km3*q_salmon), 
            cpue_var =var(JuvCatchWtCPUE_kg_km3*q_salmon),
            strata_num=n(),
            volswept=mean(SumOfOfficialVolumeSwept_km3)) %>%  
  ungroup()

## (A.2) assign the total volume (V) per stratum
mu_cpue_chinook$strata_vol <-with(mu_cpue_chinook, ifelse(STRATUM==504, 47.04, 
                                                    ifelse(STRATUM==505, 73.44, 
                                                           ifelse(STRATUM==506,24.96,
                                                                  ifelse(STRATUM==507, 40.32,
                                                                         ifelse(STRATUM==508, 37.92,
                                                                                ifelse(STRATUM==509, 44.64,
                                                                                       ifelse(STRATUM==510, 91.20, 121.92))))))))
##(B) Calculate annual stratum biomass and variance
##(B.1) Biomass
mu_cpue_chinook$strata_biomass = (mu_cpue_chinook$strata_cpue * mu_cpue_chinook$strata_vol) 

## (B.2) Variance
mu_cpue_chinook$strata_variance = mu_cpue_chinook$strata_vol*(mu_cpue_chinook$strata_vol - mu_cpue_chinook$volswept)*(mu_cpue_chinook$cpue_var/mu_cpue_chinook$strata_num)

## (C) sum up the biomass and variance estimates across all strata
biomass_chinook <-mu_cpue_chinook %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_biomass=sum(strata_biomass),
            annual_biomass_variance=sum(strata_variance),
            annual_biomass_num=sum(strata_num))

##(D) Then based on the annual biomass variance across all strata (in C above), calculate the standard deviation, CI and CV
biomass_chinook$annual_biomass_sd<-sqrt(biomass_chinook$annual_biomass_variance)
biomass_chinook$annual_biomass_cv<-biomass_chinook$annual_biomass_sd/biomass_chinook$annual_biomass
biomass_chinook$annual_biomass_se<-biomass_chinook$annual_biomass_sd/(sqrt(biomass_chinook$annual_biomass_num))
biomass_chinook$annual_biomass_LCI<-ifelse(((biomass_chinook$annual_biomass_sd*(-1.96)+biomass_chinook$annual_biomass))<0,0,(biomass_chinook$annual_biomass_sd*(-1.96)+biomass_chinook$annual_biomass)) ##if LCI<0 then assign zero
biomass_chinook$annual_biomass_UCI<-biomass_chinook$annual_biomass_sd*(1.96)+biomass_chinook$annual_biomass
biomass_chinook$species_code <- 124

## (E) Ta da!
write.table(biomass_chinook, file="biomass_estimates_fromR.v3.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",") 



#--------------------------herring -- based on q_herring and TotalCatchWtCPUE

#(A) calculate the biomass by year and day/night for each strata
## (A.1) summarize cpue and it's variance both accounting for catchability, count number of tows (n), and get the swept_volume (v) by strata
mu_cpue_herring <- cpue_herring %>%
  group_by(TRIP_YEAR, STRATUM, DayNight) %>%
  summarise(strata_cpue =mean(TotalCatchWtCPUE_kg_km3*q_herring), 
            cpue_var =var(TotalCatchWtCPUE_kg_km3*q_herring),
            strata_num=n(),
            volswept=mean(SumOfOfficialVolumeSwept_km3)) %>%  
  ungroup()

## (A.2) assign the total volume (V) per stratum
mu_cpue_herring$strata_vol <-with(mu_cpue_herring, ifelse(STRATUM==504, 47.04, 
                                                    ifelse(STRATUM==505, 73.44, 
                                                           ifelse(STRATUM==506,24.96,
                                                                  ifelse(STRATUM==507, 40.32,
                                                                         ifelse(STRATUM==508, 37.92,
                                                                                ifelse(STRATUM==509, 44.64,
                                                                                       ifelse(STRATUM==510, 91.20, 121.92))))))))
##(B) Calculate annual stratum biomass and variance
##(B.1) Biomass
mu_cpue_herring$strata_biomass = (mu_cpue_herring$strata_cpue * mu_cpue_herring$strata_vol) 

## (B.2) Variance
mu_cpue_herring$strata_variance = mu_cpue_herring$strata_vol*(mu_cpue_herring$strata_vol - mu_cpue_herring$volswept)*(mu_cpue_herring$cpue_var/mu_cpue_herring$strata_num)

## (C) sum up the biomass and variance estimates across all strata
biomass_herring <-mu_cpue_herring %>%
  group_by(TRIP_YEAR, DayNight) %>%
  summarise(annual_biomass=sum(strata_biomass),
            annual_biomass_variance=sum(strata_variance),
            annual_biomass_num=sum(strata_num))

##(D) Then based on the annual biomass variance across all strata (in C above), calculate the standard deviation, CI and CV
biomass_herring$annual_biomass_sd<-sqrt(biomass_herring$annual_biomass_variance)
biomass_herring$annual_biomass_cv<-biomass_herring$annual_biomass_sd/biomass_herring$annual_biomass
biomass_herring$annual_biomass_se<-biomass_herring$annual_biomass_sd/(sqrt(biomass_herring$annual_biomass_num))
biomass_herring$annual_biomass_LCI<-ifelse(((biomass_herring$annual_biomass_sd*(-1.96)+biomass_herring$annual_biomass))<0,0,(biomass_herring$annual_biomass_sd*(-1.96)+biomass_herring$annual_biomass)) ##if LCI<0 then assign zero
biomass_herring$annual_biomass_UCI<-biomass_herring$annual_biomass_sd*(1.96)+biomass_herring$annual_biomass
biomass_herring$species_code <- 96

## (E) Ta da!
write.table(biomass_herring, file="biomass_estimates_fromR.v3.csv", append=TRUE, row.names=FALSE, col.names=FALSE, sep=",") 




