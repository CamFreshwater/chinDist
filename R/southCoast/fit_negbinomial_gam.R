## Multinomial model fit
# January 20, 2020
# Fit negative binomial model to catch data
# Similar to fit_negbinomial but glmmTMB and gamm4 replaced with gam
# after issues with convergence

library(mgcv)
library(tidyverse)

dailyCatch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  filter(!is.na(cpue)) %>% 
         #remove areas w/ very gappy data
         #area > 100) %>% 
  mutate(area = as.factor(area),
         #condense gappy months necessary for convergence when estimating
         #month:area interactions or random month slopes by area 
         month = as.factor(month),
         month_n = as.numeric(month),
         year = as.factor(year),
         z_eff = as.numeric(scale(boatDays)),
         z_eff2 = z_eff^2) %>% 
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
  geom_point(aes(x = z_eff2, y = catch)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~area, scales = "free")

table(dailyCatch$month, dailyCatch$reg)


## Fit model -------------------------------------------------------------------
# Goal is to make stat area predictions by month, while accounting for effort
# and variation among years. Consider making predictions for NW and SWVI, using
# stat area as a random effect only

# Compare various models
# Model 1 = month intercepts w/ RE intercepts for year and area
nb <- gam(catch ~ s(eff) + month + s(year, bs = "re") + s(area, bs = "re"), 
            data = dailyCatch,
            family = nb)
# Adds a cyclic crs to month
nb_cc <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + s(year, bs = "re") + 
                 s(area, bs = "re"), 
            data = dailyCatch,
            family = nb,
            knots = list(month_n = c(1, 12)))
# Adds a tensor product between month and area
nb_cc_tp <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + 
                    te(month_n, area, bs = c("cc", "re"), m = 2) + 
                    s(year, bs = "re"), 
               data = dailyCatch,
               family = nb,
               knots = list(month_n = c(1, 12)))
# Uses poisson instead of neg binom
pois_cc_tp <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + 
                    te(month_n, area, bs = c("cc", "re"), m = 2) + 
                    s(year, bs = "re"), 
                  data = dailyCatch,
                  family = poisson,
                  knots = list(month_n = c(1, 12)))
# Adds a tensor product between effort and month
nb_cc_tp2 <- gam(catch ~ te(z_eff, month_n, bs = c("tp", "cc"), m = 2) + 
                    t2(month_n, area, bs = c("cc", "re"),
                             m = 2, full = TRUE) + 
                    s(year, bs = "re"), 
                  data = dailyCatch,
                  family = nb,
                  knots = list(month_n = c(1, 12)))

AIC(nb); AIC(nb_cc); AIC(nb_cc_tp); AIC(pois_cc_tp); AIC(nb_cc_tp2)
# poisson favored but strong evidence of overdisp so stick w/ neg binomial

gam.vcomp(nb_cc_tp)


## Check for overdispersion
# root_pois <- countreg::rootogram(pois_cc_tp, style = "hanging", plot = TRUE)
# root_nb   <- countreg::rootogram(nb_cc_tp, style = "hanging", plot = TRUE)
# 
# sum(residuals(pois_cc_tp, type = "pearson")^2) / df.residual(pois_cc_tp)
# sum(residuals(nb_cc_tp, type = "pearson")^2) / df.residual(nb_cc_tp)
# definitely necessary to use negative binomial



## Check model -----------------------------------------------------------------

plot_fits(nb_cc_tp, dat = dailyCatch, exclude = FALSE, y = "resid")
plot_fits(nb_cc_tp, dat = dailyCatch, exclude = TRUE, y = "obs")

plot_fits <- function(mod, dat, exclude = TRUE, y = "resid", nbin = TRUE) {
  fit_dat <- broom.mixed::augment(mod, data = dat)
  
  if (exclude == FALSE) {
    if (y == "resid") {
      q <- ggplot(fit_dat, aes(.fitted, .resid))  
    }
    if (y == "obs") {
      q <- ggplot(fit_dat, aes(.fitted, catch))  
    }
    p <- q + geom_point(aes(colour = area)) +
      # geom_smooth() +
      ggsidekick::theme_sleek() +
      facet_wrap(~month, scales = "free")
  }
  
  if (exclude == TRUE) {
    yy <- unique(dat$year)[1]
    aa <- unique(dat$area)[1]
    
    newd <- transform(dat, year = yy)
                      #, area = aa)
    pred <- predict(mod, newd, exclude = c("s(year)"))
    #transform predictions if using nbinom
    if(nbin == TRUE) {
      pred <- exp(pred)
    }
    fit_dat$.fitted0 <- pred
    fit_dat$.resid0 <- fit_dat$catch - fit_dat$.fitted0
    
    if (y == "resid") {
      q <- ggplot(fit_dat, aes(.fitted0, .resid0))  
    }
    if (y == "obs") {
      q <- ggplot(fit_dat, aes(.fitted0, catch))  
    }
    p <- q + 
      geom_point(aes(colour = area)) +
      # geom_smooth() +
      ggsidekick::theme_sleek() +
      facet_wrap(~month, scales = "free")
  }
  return(p)
}



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




