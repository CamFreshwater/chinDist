## Multinomial model fit
# January 20, 2020
# Fit negative binomial model to catch data

library(glmmTMB)
library(tidyverse)

weekCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                                 "weeklyCatch_WCVI.csv"), 
                      stringsAsFactors = TRUE) %>% 
  mutate(
    statArea = as.factor(statArea),
    month = as.factor(month),
    year = as.factor(year),
    weeklyCPUE = weeklyCatch / weeklyEffort)

monthC <- weekCatch %>%
  group_by(month, year, statArea) %>% 
  summarize(catch = sum(weeklyCatch),
            eff = sum(weeklyEffort),
            cpue = catch / eff) %>% 
  filter(!is.na(cpue),
         #remove inshore regions with highly variable catches
         #statArea %in% c("24", "25", "26"),
         !cpue > 150) %>% #remove temp outlier 
  ungroup() %>% 
  mutate(
    reg = case_when(
      statArea %in% c("123", "124", "23", "24") ~ "SWVI",
      TRUE ~ "NWVI"
    ),
    reg = as.factor(reg)
  )


## Visualize
ggplot(monthC, aes(x = month, y = cpue)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  facet_wrap(~statArea)

ggplot(monthC) +
  geom_point(aes(x = catch, y = eff)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~month, scales = "free")


## Fit model -------------------------------------------------------------------
# Goal is to make stat area predictions by month, while accounting for effort
# and variation among years. Consider making predictions for NW and SWVI, using
# stat area as a random effect only

mod1 <- glmmTMB(catch ~ eff + month * reg + (1|statArea) + (1|year), 
                family = nbinom2,
                data = monthC)


gmod_gA_L_NB2 <- glmmadmb(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                          family="nbinom",data=Gdat)
gmod_gA_L_NB1 <- glmmadmb(shells~prev+offset(log(Area))+factor(year)+(1|Site),
                          family="nbinom1",data=Gdat)