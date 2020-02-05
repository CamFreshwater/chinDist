## Multinomial model fit
# January 20, 2020
# Fit negative binomial model to catch data

library(glmmTMB)
library(tidyverse)

dailyCatch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.RDS"))
monthC <- dailyCatch %>% 
  group_by(month, year, area, catchReg) %>% 
  summarize(catch = sum(catch),
            eff = sum(boatDays),
            cpue = catch / eff) %>% 
  filter(!is.na(cpue),
         #remove inshore regions with highly variable catches
         #statArea %in% c("24", "25", "26"),
         !cpue > 150) %>% #remove temp outlier 
  ungroup() %>% 
  rename(reg = catchReg) %>% 
  mutate(month = as.factor(month),
         year = as.factor(year),
         z_eff = scale(eff))


## Visualize
ggplot(monthC, aes(x = month, y = cpue)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  facet_wrap(~area)

ggplot(monthC) +
  geom_point(aes(x = catch, y = eff)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~month, scales = "free")


## Fit model -------------------------------------------------------------------
# Goal is to make stat area predictions by month, while accounting for effort
# and variation among years. Consider making predictions for NW and SWVI, using
# stat area as a random effect only

#nbinom2 supported rel to 1 based on aic scores
mod2 <- glmmTMB(catch ~ z_eff + month + reg + (1|area) + (1|year), 
                family = nbinom2,
                data = monthC)



gmod_gA_L_NB2 <- glmmadmb(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                          family="nbinom",data=Gdat)
gmod_gA_L_NB1 <- glmmadmb(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                          family="nbinom1",data=Gdat)