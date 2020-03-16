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
  mutate(catchReg = as.factor(catchReg),
         area = as.factor(area),
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
# nb <- gam(catch ~ s(eff) + month + s(year, bs = "re") + s(area, bs = "re"), 
#             data = dailyCatch,
#             family = nb)
# # Adds a cyclic crs to month
# nb_cc <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + 
#                #s(year, bs = "re") + 
#                  s(area, year, bs = "re"), 
#             data = dailyCatch,
#             family = nb,
#             knots = list(month_n = c(1, 12)))
# Adds a tensor product between month and area
nb_cc_tp <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") +
                    te(month_n, area, bs = c("cc", "re"), m = 2) +
                    s(area, year, bs = "re"),
               data = dailyCatch,
               family = nb,
               knots = list(month_n = c(1, 12)))
# Uses poisson instead of neg binom
# pois_cc_tp <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + 
#                     te(month_n, area, bs = c("cc", "re"), m = 2) + 
#                     s(year, bs = "re"), 
#                   data = dailyCatch,
#                   family = poisson,
#                   knots = list(month_n = c(1, 12)))
# Adds a tensor product between month and area, then nests area within year
nb_cc_tp2 <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + 
                   s(month_n, area, bs = "fs", xt = list(bs = "cc")) +
                   s(area, year, bs = "re"), 
                 data = dailyCatch,
                 family = nb,
                 knots = list(month_n = c(1, 12)))

# Adds a tensor product between month, area and year
nb_cc_tp3 <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + 
                  t2(month_n, area, year, bs = c("cc", "re", "re"), m = 2,
                     full = TRUE), 
                data = dailyCatch,
                family = nb,
                knots = list(month_n = c(1, 12)))

# Adds a tensor product between month/area and month/area/year
nb_cc_tp4 <- gam(catch ~ s(z_eff) + s(month_n, bs = "cc") + 
                   te(month_n, area, bs = c("cc", "re"), m = 2) +
                   t2(month_n, area, year, bs = c("cc", "re", "re"), m = 2,
                      full = TRUE), 
                 data = dailyCatch,
                 family = nb,
                 knots = list(month_n = c(1, 12)))

# More complex model to estimate fixed regional effect (consider using, but not
# for now)
nb_cc_tp5 <- gam(catch ~ s(z_eff) + #te(month_n, bs = "cc") + 
                   s(month_n, reg, bs = "fs", xt = list(bs = "cc")) +
                   t2(month_n, area, year, bs = c("cc", "re", "re"), m = 2,
                      full = TRUE), 
                 data = dailyCatch,
                 family = nb,
                 knots = list(month_n = c(1, 12)))

# nb_cc_tp3_fix <- gam(catch ~ s(z_eff, k = 20) + s(month_n, bs = "cc", k = 10) + 
#                    t2(month_n, area, year, bs = c("cc", "re", "re"), k = 10, m = 2,
#                       full = TRUE), 
#                  data = dailyCatch,
#                  family = nb,
#                  knots = list(month_n = c(1, 12)))




# AIC(nb); AIC(nb_cc); 
AIC(nb_cc_tp); AIC(nb_cc_tp2); AIC(nb_cc_tp3); AIC(nb_cc_tp4); AIC(nb_cc_tp5)
# AIC(nb_cc_tp3_fix)
# poisson favored but strong evidence of overdisp so stick w/ neg binomial

gam.check(nb_cc_tp4)
gam.check(nb_cc_tp3)

k.check(nb_cc_tp2)
k.check(nb_cc_tp3_fix)
# some evidence that default k is insufficient

gam.vcomp(nb_cc_tp3)


## Check for overdispersion
# root_pois <- countreg::rootogram(pois_cc_tp, style = "hanging", plot = TRUE)
# root_nb   <- countreg::rootogram(nb_cc_tp, style = "hanging", plot = TRUE)
# 
# sum(residuals(pois_cc_tp, type = "pearson")^2) / df.residual(pois_cc_tp)
# sum(residuals(nb_cc_tp, type = "pearson")^2) / df.residual(nb_cc_tp)
# definitely necessary to use negative binomial



## Check model -----------------------------------------------------------------

plot_fits <- function(mod, dat, exclude = TRUE, y = "resid") {
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
    pred <- predict(mod, newd, exclude = c("s(year)"))
    #transform predictions if using nbinom
    fit_dat$.fitted0 <- exp(pred)
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

plot_fits(nb_cc_tp3, dat = dailyCatch, exclude = FALSE, y = "resid")
plot_fits(nb_cc_tp3, dat = dailyCatch, exclude = TRUE, y = "obs")


## Generate predictions with z_eff zeroed
newDat <- expand.grid(area = unique(dailyCatch$area), 
                      month_n = unique(dailyCatch$month_n), 
                      year =  unique(dailyCatch$year)) %>% 
  mutate(z_eff = 0) 
# only include if estimating fixed effects for reg %>% 
  # left_join(., 
  #           dailyCatch %>% 
  #             select(reg, area) %>% 
  #             distinct(),
  #           by = "area")


# Make various predictions
# head(predict(nb_cc_tp, newdata = newDat)) 
# head(predict(nb_cc_tp, newdata = newDat, exclude = "s(area,year)")) 
# head(predict(nb_cc_tp2, newdata = newDat)) 
# head(predict(nb_cc_tp2, newdata = newDat, exclude = "s(area,year)")) 
head(predict(nb_cc_tp3, newdata = newDat)) 
tt <- predict.gam(nb_cc_tp3, newdata = newDat, se.fit = T) %>% 
  bind_cols()

pred_fe <- predict(nb_cc_tp3, newdata = newDat, 
                   exclude = "t2(month_n,area,year)")


# Plot different models
newDat %>% 
  mutate(tp1 = as.numeric(predict(nb_cc_tp, newdata = newDat)),
         tp2 = as.numeric(predict(nb_cc_tp2, newdata = newDat)),
         tp3 = as.numeric(predict(nb_cc_tp3, newdata = newDat))) %>%
  pivot_longer(cols = tp1:tp3,
               names_to = "model",
               values_to = "preds") %>% 
  ggplot(.) +
  geom_boxplot(aes(x = as.factor(month_n), y = preds, fill = model)) +
  ggsidekick::theme_sleek()
#relatively similar median estimates

# compare two most saturated models with random effects included
newDat %>% 
  mutate(preds = as.numeric(predict(nb_cc_tp4, newdata = newDat))) %>%
  ggplot(.) +
  geom_point(aes(x = as.factor(month_n), y = preds)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~area)
# newDat %>% 
#   mutate(preds = as.numeric(predict(nb_cc_tp5, newdata = newDat))) %>%
#   ggplot(.) +
#   geom_point(aes(x = as.factor(month_n), y = preds)) +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~area)


# Extract fixed effects estimates from top model
# Plot pooling across areas/years
preds <- predict.gam(nb_cc_tp4, newdata = newDat, 
                     exclude = c("t2(month_n,area,year)"),
                     se = T)

mu_dat1 <- newDat %>% 
  mutate(mu.pred = as.numeric(preds$fit),
         se.pred = as.numeric(preds$se.fit),
         low = mu.pred - (qnorm(0.975) * se.pred),
         high = mu.pred + (qnorm(0.975) * se.pred)) %>% 
  select(-year) %>% 
  distinct()

ggplot(mu_dat1, aes(x = as.factor(month_n), y = mu.pred)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~area)


