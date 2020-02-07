## Multinomial model fit
# January 20, 2020
# Fit negative binomial model to catch data

library(glmmTMB)
library(tidyverse)

dailyCatch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.RDS")) %>% 
  filter(!is.na(sumCPUE)) %>% 
  mutate(area = as.factor(area),
         month = case_when(
           month %in% c("6", "7") ~ "6_7",
           TRUE ~ as.character(month)
         ),
         month = as.factor(month),
         year = as.factor(year),
         z_eff = scale(boatDays)) %>% 
  rename(eff = boatDays,
         cpue = sumCPUE,
         reg = catchReg)

## Visualize
ggplot(dailyCatch, aes(x = month, y = cpue)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  facet_wrap(~reg)

ggplot(dailyCatch, aes(x = month, y = eff)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  facet_wrap(~reg)

ggplot(dailyCatch) +
  geom_point(aes(x = catch, y = eff)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~month, scales = "free")

table(dailyCatch$month, dailyCatch$reg)

## Fit model -------------------------------------------------------------------
# Goal is to make stat area predictions by month, while accounting for effort
# and variation among years. Consider making predictions for NW and SWVI, using
# stat area as a random effect only

#nbinom2 supported rel to 1 based on aic scores
mod2 <- glmmTMB(cpue ~ month + reg + (1|area) + (1|year), 
                    family = nbinom2,
                    data = dailyCatch)
mod2_int <- glmmTMB(cpue ~ (month:reg) + (1|area) + (1|year), 
                family = nbinom2,
                data = dailyCatch)



library(lme4)
mod3 <- glmer.nb(catch ~ z_eff + month + reg + (1|area) + (1|year), 
         data = monthC)


mod <- lmer(FL ~ zDOY + SEX + AGE + (1 | YEAR), data = FullDataset,
            REML = FALSE)



gmod_gA_L_NB2 <- glmmadmb(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                          family="nbinom",data=Gdat)
gmod_gA_L_NB1 <- glmmadmb(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                          family="nbinom1",data=Gdat)