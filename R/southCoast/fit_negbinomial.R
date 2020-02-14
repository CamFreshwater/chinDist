## Multinomial model fit
# January 20, 2020
# Fit negative binomial model to catch data

library(lme4)
library(glmmTMB)
library(tidyverse)

dailyCatch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  filter(!is.na(cpue),
         #remove areas w/ very gappy data
         area > 100) %>% 
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
  geom_point(aes(x = z_eff^2, y = catch)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~area, scales = "free")

table(dailyCatch$month, dailyCatch$reg)

## Fit model -------------------------------------------------------------------
# Goal is to make stat area predictions by month, while accounting for effort
# and variation among years. Consider making predictions for NW and SWVI, using
# stat area as a random effect only


# converges but SEs can't be estimated
# mod2 <- glmmTMB(catch ~ z_eff + month:reg + (1|area) + (1|year),
#                 family = nbinom2,
#                 data = dailyCatch)

# fit regions separately since different data available and interaction model 
# can't be fit
# negative binomial has very poor fits due to non-linear relationship
dc_north <- dailyCatch %>% 
  filter(reg == "NWVI") %>% 
  mutate(z_eff = scale(eff))
mod2_n <- glmmTMB(catch ~ (z_eff)^2 + month + (1|area) + (1|year),
                  family = nbinom2, data = dc_north)

no_zero_n <- dc_north %>% 
  filter(!cpue == "0")
mod2_n_g <- glmmTMB(cpue ~ month + (1|area) + (1|year), 
                    family = Gamma(link = "log"),
                    data = no_zero_n)

dc_south <- dailyCatch %>% 
  filter(reg == "SWVI") %>% 
  mutate(z_eff = scale(eff))
mod2_s <- glmmTMB(catch ~ (z_eff^2) + month + (1|area) + (1|year),
                  family = nbinom2, data = dc_south)

#glmer gives similar results, but much slower and has warnings
mod2b <- glmer.nb(catch ~ z_eff + month + (1|area) + (1|year), 
                        #family = nbinom2,
                        data = dailyCatch)


## Experiment with gamms given evidence of non-linearities between effort and 
# catch which result in skewed residuals
library(gamm4)
mod3_n <- gamm4(catch ~ s(z_eff) + month, data=dc_north, 
              random = ~(1|area) + (1|year),
              family = nbinom2)
mod3_n_nest <- gamm4(catch ~ s(z_eff) + month, data = dc_north, 
                    random = ~(1|year:area),
                    family = nbinom2)

no_zero_n <- dc_north %>% 
  filter(!month == "11")
mod3_y <- gamm4(catch ~ month, data = no_zero_n, 
                     random = ~(1|year),
                     family = nbinom2)
summary(mod3_y$mer)
mod3_a <- gamm4(catch ~ s(z_eff) + month, data = no_zero_n, 
                random = ~(1|area),
                family = nbinom2)



mod3_s <- gamm4(catch ~ s(z_eff) + month, data = dc_south, 
                random = ~(1|area) + (1|year),
                family = nbinom2)
mod3_s_int <- gamm4(catch ~ s(z_eff) + month, data = dc_south, 
                random = ~(month|area) + (1|year),
                family = nbinom2)
plot(mod3_n$gam)
plot(mod3_n$mer)
plot(mod3_n_int$gam)
plot(mod3_n_int$mer)

summary(mod3_n_nest$gam)


set.seed(0)
dat <- gamSim(1,n=400,scale=2) ## simulate 4 term additive truth
## Now add 20 level random effect `fac'...
dat$fac <- fac <- as.factor(sample(1:20,400,replace=TRUE))
dat$y <- dat$y + model.matrix(~fac-1)%*%rnorm(20)*.5

br <- gamm4(y~s(x0)+x1+s(x2),data=dat,random=~(1|fac))
summary(br$mer)

## Check model -----------------------------------------------------------------

# Check model residuals w/ DHARMa 
mod2s_resid <- DHARMa::simulateResiduals(mod3_n$mer)
plot(mod2s_resid)
#some issues w/ residuals vs. predicted


aa_gamm <- broom.mixed::augment(mod3_n$mer, data = dc_north)
aa_nb <- broom.mixed::augment(mod2_n, data = dc_north)
aa_g <- broom.mixed::augment(mod2_n_g, data = no_zero_n)
# aa_s <- broom.mixed::augment(mod3_s$mer, data = dc_south)
ggplot(aa_nb, aes(.fitted, .resid)) +
  geom_line(aes(group = area), colour = "gray") +
  # geom_line(aes(group = year), colour = "gray") +
  geom_point(aes(colour = month)) +
  geom_smooth() +
  ggsidekick::theme_sleek()


aa_gamm$.fitted0 <- predict(mod3_n$gam, newdata = transform(dc_north, year = NA, 
                                                   area = NA), 
                       type = "response")
aa_gamm$.resid0 <- dc_north$catch - aa_gamm$.fitted0
aa_g$.fitted0 <- predict(mod2_n_g,
                         newdata = transform(no_zero_n, year = NA, area = NA), 
                       type = "response")
aa_g$.resid0 <- no_zero_n$catch - aa_g$.fitted0
aa_nb$.fitted0 <- predict(mod2_n,
                         newdata = transform(dc_north, year = NA, area = NA), 
                         type = "response")
aa_nb$.resid0 <- dc_north$catch - aa_nb$.fitted0

ggplot(aa_nb, aes(.fitted0,.resid0)) + 
  geom_line(aes(group=year), colour="gray") + 
  geom_point(aes(colour=month)) + 
  geom_smooth()



# Check predictions of catch gamm's with neg. binomial vs. cpue glmm w/ gamma
# glmm
new_dat <- data.frame(month = levels(dc_north$month),
                      z_eff = 0)
xx <- model.matrix(lme4::nobars(formula(mod2_n_g)[-2]), new_dat)
cond_betas <- fixef(mod2_n_g)$cond
preds <- xx %*% cond_betas

# gamm
newdat_gamm <-  data.frame(month = levels(dc_north$month),
                           z_eff = 0)
preds <- predict(mod3_n$gam, newdata = newdat_gamm, se.fit = TRUE) %>% 
  cbind(newdat_gamm, .) %>% 
  mutate(low = fit - (1.96 * se.fit),
         up = fit + (1.96 * se.fit),
         month = fct_relevel(as.factor(month), "10", after = 10, "11", 
                             after = 11, "12", after = 12))

ggplot(preds, aes(x = month, y = fit)) +
  geom_pointrange(aes(ymin = low, ymax = up), shape = 21,
                  position = position_dodge(0.9)) +
  # facet_wrap(~year) +
  ggsidekick::theme_sleek()

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




