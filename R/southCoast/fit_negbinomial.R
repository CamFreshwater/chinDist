## Multinomial model fit
# January 20, 2020
# Fit negative binomial model to catch data

library(lme4)
library(glmmTMB)
library(tidyverse)

dailyCatch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  filter(!is.na(cpue)) %>% 
  mutate(area = as.factor(area),
         #condense gappy months necessary for convergence when estimating
         #month:area interactions or random month slopes by area 
         month_c = case_when(
           month %in% c("6", "7") ~ "6_7",
           TRUE ~ as.character(month)
         ),
         month_c = as.factor(month_c),
         month = as.factor(month),
         year = as.factor(year),
         z_eff = as.numeric(scale(boatDays))) %>% 
  rename(eff = boatDays,
         reg = catchReg)

## Visualize
ggplot(dailyCatch, aes(x = month, y = cpue)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  facet_wrap(~reg)

ggplot(dailyCatch, aes(x = month, y = eff)) +
  geom_boxplot() +
  ggsidekick::theme_sleek() +
  facet_wrap(~area)

ggplot(dailyCatch) +
  geom_point(aes(x = catch, y = eff)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~month, scales = "free")

table(dailyCatch$month, dailyCatch$reg)

## Fit model -------------------------------------------------------------------
# Goal is to make stat area predictions by month, while accounting for effort
# and variation among years. Consider making predictions for NW and SWVI, using
# stat area as a random effect only


#nbinom1 does not converge
# mod1_int <- glmmTMB(catch ~ z_eff + (month:reg) + (1|area) + (1|year), 
#                     family = nbinom1,
#                     data = dailyCatch)

# mod2 <- glmmTMB(catch ~ z_eff + month + (1|area) + (1|year),
#                 family = nbinom2,
#                 data = dailyCatch)
# fit separately since different data available and interaction model can't be fit
mod2_n <- glmmTMB(catch ~ z_eff + month + (1|area) + (1|year),
                  family = nbinom2,
                  data = dailyCatch %>% filter(reg == "NWVI"))
mod2_s <- glmmTMB(catch ~ z_eff + month + (1|area) + (1|year),
                  family = nbinom2,
                  data = dailyCatch %>% filter(reg == "SWVI"))

#glmer gives similar results, but much slower and has warnings
mod2b <- glmer.nb(catch ~ z_eff + month + (1|area) + (1|year), 
                        #family = nbinom2,
                        data = dailyCatch)

# fails to converge due to data gaps
# mod2_rs <- glmmTMB(catch ~ z_eff + month + (1 + month|area) + (1|year), 
#                 family = nbinom2,
#                 data = dailyCatch)


mod2_reg <- glmmTMB(catch ~ z_eff + month + reg + (1|area) + (1|year),
                    family = nbinom2,
                    data = dailyCatch)

mod3_nest <- glmmTMB(catch ~ z_eff + (month:reg) + (1|reg:area) + (1|year), 
                    family = nbinom2,
                    data = dailyCatch)

bbmle::AICtab(mod2, mod2_int, mod3_nest)


## Inference -------------------------------------------------------------------

# Check model residuals
mod2_resid <- DHARMa::simulateResiduals(mod2)
plot(mod2_resid)


# Check predictions
new_dat <- data.frame(month = levels(dailyCatch$month),
                      z_eff = 0)
xx <- model.matrix(lme4::nobars(formula(mod2)[-2]), new_dat)
cond_betas <- fixef(mod2)$cond
preds <- xx %*% cond_betas


# Simulate from model
sims <- simulate(mod2_s, seed = 1, nsim = 1000)
simdat <- map(sims, function(count) {
  cbind(count, 
        dailyCatch %>% 
          filter(reg == "SWVI") %>% 
          select(month, area, year)) %>%
    group_by(month) %>% 
    summarize(mu = mean(count))
}) %>% 
  bind_rows() 



simdatlist=lapply(sims, function(count){
  cbind(count, Salamanders[,c('site', 'mined', 'spp')])
})
simdatsums=lapply(simdatlist, function(x){
  ddply(x, ~spp+mined, summarize,
        absence=mean(count==0),
        mu=mean(count))
})
ssd=do.call(rbind, simdatsums)



zinbm3 = glmmTMB(count~spp * mined +(1|site), zi=~spp * mined, Salamanders, 
                 family=nbinom2)




