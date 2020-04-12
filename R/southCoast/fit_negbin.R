## Neg bin model fit
# April 6, 2020
# Fit neg binomial model to catch data from WCVI troll fishery
# Persistently underperforms compared to CPUE model w/ tweedie; don't use for 
# now; should probably include a zero inflated or heavy-tailed version

library(tidyverse)
library(TMB)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  filter(!is.na(boatDays),
         !boatDays == "0",
         !month == "7") %>% 
  rename(eff = boatDays) %>% 
  mutate(reg = factor(catchReg),
         area = as.factor(area),
         month = as.factor(month),
         month_n = as.numeric(month),
         year = as.factor(year),
         log_cpue = log(cpue + 0.0001),
         eff_z = as.numeric(scale(eff)),
         eff_z2 = eff_z^2,
         eff_z3 = eff_z^3) %>% 
  arrange(catchReg, month)
hist(catch$cpue)
hist(catch$catch)
hist(catch$eff_z2)
plot(catch ~ eff, data = catch)

# Prep data to pass to model
yr_vec <- as.numeric(as.factor(as.character(catch$year))) - 1

# model matrix for fixed effects
fix_mm <- model.matrix(~ -1 + (catchReg : month) + eff_z + eff_z2, catch)
# model matrix for predictions
# mm_pred1 <- fix_mm[,-c(1,2)] %>% 
#   unique() 
# mm_pred <- cbind(mm_pred1,
#                  eff_z = rep(0, n = nrow(mm_pred1)),
#                  eff_z2 = rep(0, n = nrow(mm_pred1)))

# Factor key of unique combinations to generate predictions
fac_key <- catch %>%
  select(catchReg, month) %>%
  distinct() %>%
  mutate(facs = paste(catchReg, month, sep = "_"),
         # facs = fct_relevel(facs, c(as.character(facs))),
         facs = fct_reorder2(facs, catchReg, 
                             desc(as.numeric(as.character(month)))),
         facs_n = as.numeric(as.factor(facs))
  ) %>%
  arrange(facs_n)
mm_pred <- model.matrix(~ facs - 1, fac_key) %>% 
  cbind(eff_z = rep(0, n = nrow(.)),
        eff_z2 = rep(0, n = nrow(.)),
        .)

data <- list(y1_i = catch$catch,
             X1_ij = fix_mm,
             factor1k_i = yr_vec,
             nk1 = length(unique(yr_vec)),
             X1_pred_ij = mm_pred
)

# Fit simple model to initialize tmb
# m1 <- lm(log(catch + 0.0001) ~ catchReg + month + eff_z + eff_z2, data = catch)
m1 <- lm(log(catch + 0.0001) ~ -1 + (catchReg : month) + eff_z + eff_z2, 
         data = catch)
m2 <- glmmTMB::glmmTMB(catch ~ catchReg + month + eff_z, data = catch, 
                       family = glmmTMB::nbinom2(link = "log")
)
m2a <- glmmTMB::glmmTMB(catch ~ catchReg + month + eff_z + eff_z2, data = catch, 
                       family = glmmTMB::nbinom2(link = "log")
)
m2b <- glmmTMB::glmmTMB(catch ~ -1 + (catchReg : month) + eff_z + eff_z2, data = catch, 
                        family = glmmTMB::nbinom2(link = "log")
)
summary(m2)
summary(m2a)
summary(m2b)

parameters = list(
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
  log_phi = log(1.5),
  z1_k = rep(0, length(unique(yr_vec))),
  log_sigma_zk1 = log(0.25)
)


# FIT --------------------------------------------------------------------------
# Compile
compile(here::here("R", "southCoast", "tmb", "negbin_1re.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", "negbin_1re")))

obj <- MakeADFun(data, parameters, random = c("z1_k"), 
                 DLL = "negbin_1re")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr

saveRDS(ssdr, here::here("generatedData", "model_fits", "negbin_ssdr.RDS"))


# PREDICTIONS ------------------------------------------------------------------

ssdr <- readRDS(here::here("generatedData", "model_fits", "negbin_ssdr.RDS"))

## Plot predictions
log_pred_fe <- ssdr[rownames(ssdr) %in% "log_prediction", ]
pred_ci <- data.frame(log_pred_est = log_pred_fe[ , "Estimate"],
                      log_pred_se =  log_pred_fe[ , "Std. Error"]) %>%
  mutate(facs_n = fac_key$facs_n) %>% 
  mutate(log_pred_low = log_pred_est + (qnorm(0.025) * log_pred_se),
         log_pred_up = log_pred_est + (qnorm(0.975) * log_pred_se),
         pred_est = exp(log_pred_est),
         pred_se = exp(log_pred_se),
         pred_low = exp(log_pred_low),
         pred_up = exp(log_pred_up)) %>%
  left_join(., fac_key, by = "facs_n") 

ggplot() +
  geom_boxplot(data = catch, 
               aes(x = as.factor(month), y = log(catch + 0.001)), alpha = 0.2) +
  geom_pointrange(data = pred_ci, aes(x = as.factor(month), y = log_pred_est,
                                      ymin = log_pred_low, 
                                      ymax = log_pred_up)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~catchReg)

pred_catch <- ggplot() +
  geom_boxplot(data = catch, aes(x = as.factor(month), y = catch / eff_z), 
               alpha = 0.2) +
  geom_pointrange(data = pred_ci, aes(x = as.factor(month), y = pred_est,
                                      ymin = pred_low, ymax = pred_up)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~catchReg)

ggplot(catch) +
  geom_point(aes(x = eff_z, y = catch), 
             alpha = 0.2) +
  lims(x = c(0, 4)) +
  # stat_function(fun = function(x) 204 + 412*x)
  stat_function(fun = function(x) exp(4.13 + 1.873*x - 0.2662*x^2)) +
  stat_function(fun = function(x) exp(3.99 + 1.275*x), color = "blue") +
  stat_function(fun = function(x) exp(4.34 + 1.876*x - 0.263*x^2), color = "red")


xx <- rnbinom(10000, mu = pred_ci$pred_est, size = exp(0.154))


pdf(here::here("figs", "model_pred", "nb_catch_pred.pdf"),
    height = 4.5, width = 6)
pred_catch
pred_catch +
  lims(y = c(0, 1000))
dev.off()
