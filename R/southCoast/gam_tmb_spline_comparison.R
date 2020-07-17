## Compare GAM and TMB spline estimates
# July 17, 2020 
# Simple fixed-effects only model without regional aggregation that compares 
# predictions from GAM and TMB spline models

library(tidyverse)
library(TMB)
library(mgcv)

# clean function 
clean_catch <- function(dat) {
  dat %>% 
    filter(!is.na(eff),
           !eff == "0") %>% 
    mutate(region = factor(region),
           area_n = as.numeric(area),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff)),
           eff_z2 = eff_z^2) %>% 
    arrange(region, month) 
}

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  rename(eff = boatDays, region = catchReg) %>% 
  clean_catch(.) %>% 
  # drop inside areas where seasonal catches not available
  filter(!area_n < 100) %>% 
  droplevels()

#recreational catch data
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "monthlyCatch_rec.RDS")) %>% 
  #group by subarea to get rid of adipose and released legal categories
  group_by(month, month_n, year, area, subarea, region, legal) %>% 
  mutate(subarea_catch = sum(mu_catch),
         subarea_eff = mean(mu_boat_trips)) %>% 
  #group by area to be consistent with commercial data
  group_by(month, month_n, year, area, region, legal) %>% 
  summarize(catch = sum(subarea_catch),
            eff = sum(subarea_eff),
            cpue = catch / eff) %>% 
  ungroup() %>% 
  mutate(temp_strata = paste(month_n, region, sep = "_"),
         region_c = as.character(region),
         region = abbreviate(region, minlength = 4))  %>% 
  clean_catch(.) %>% 
  # drop months with minimal catch estimates
  filter(!month_n < 5,
         !month_n > 9) %>% 
  droplevels()


# PREP DATA --------------------------------------------------------------------

#subset of data
catch_dat <- comm_catch %>% 
  filter(region == "NWVI",
         !month_n < 5,
         !month_n > 9)

yr_vec <- as.numeric(as.factor(as.character(catch_dat$year))) - 1

#generate model matrix based on GAM
months <- unique(catch_dat$month_n)
n_months <- length(months)
m1 <- gam(catch ~ s(month_n, bs = "tp", k = kk, by = area) +
                 eff_z + eff_z2, 
               knots = list(month_n = c(min(months), max(months))),
               data = catch_dat,
               family = nb)
fix_mm <- predict(m1, type = "lpmatrix")

# make predictive model matrix including null values for effort
pred_dat <- expand.grid(
  month_n = seq(min(catch_dat$month_n), 
                max(catch_dat$month_n),
                length.out = 50),
  area = unique(catch_dat$area),
  eff_z = 0,
  eff_z2 = 0
)
pred_mm_catch <- predict(m1, pred_dat, type = "lpmatrix")

data <- list(y1_i = catch_dat$catch,
             X1_ij = fix_mm,
             factor1k_i = yr_vec,
             nk1 = length(unique(yr_vec)),
             X1_pred_ij = pred_mm_catch
             )

parameters <- list(
  #abundance parameters
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
  log_phi = log(1.5)
)


# FIT MODEL --------------------------------------------------------------------

compile(here::here("src", "negbin_fe.cpp"))
dyn.load(dynlib(here::here("src", "negbin_fe")))


obj <- MakeADFun(data, parameters, DLL = "negbin_fe")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
ssdr <- summary(sdr)


#Predictions using MGCV
pred_y <- predict.gam(m1, pred_dat, type = "response") 
pred_catch <- pred_dat %>% 
  mutate(y = pred_y)

# TMB predictions
abund_pred <- ssdr[rownames(ssdr) %in% "pred_abund", ] 
pred_ci <- data.frame(pred_est = abund_pred[ , "Estimate"],
                      pred_se =  abund_pred[ , "Std. Error"]) %>%
  cbind(pred_dat, .) %>% 
  mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
         pred_up = pred_est + (qnorm(0.975) * pred_se)) 


ggplot(data = pred_ci, aes(x = month_n)) +
  geom_line(aes(y = pred_est#, colour = region
  )) +
  geom_ribbon(aes(ymin = pred_low, ymax = pred_up
                  #, fill = region
  ), alpha = 0.5)  +
  geom_line(data = pred_catch, aes(x = month_n, y = y), col = "red") +
  # geom_point(data = catch_dat, aes(x = month_n, y = cpue)) +
  facet_wrap(~area)
