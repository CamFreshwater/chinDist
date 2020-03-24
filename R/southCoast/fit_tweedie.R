## Tweedie model fit
# March 20, 2020
# Fit gamma model to CPUE data from WCVI troll fishery

library(tidyverse)
library(TMB)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  filter(!is.na(cpue)) %>% 
  mutate(reg = factor(catchReg),
         area = as.factor(area),
         month = as.factor(month),
         month_n = as.numeric(month),
         year = as.factor(year),
         log_cpue = log(cpue + 0.0001)) %>% 
  rename(eff = boatDays) %>% 
  arrange(catchReg, month)
hist(catch$cpue)
hist(catch$log_cpue)


# helper function for inverse logit of p to check par estimates
inv_logit <- function(x) {
  1.0 / (1.0 + exp(-x)) + 1.0
}


# Prep data to pass to model
y_obs <- catch$cpue
yr_vec <- as.numeric(as.factor(as.character(catch$year))) - 1

# model matrix for fixed effects
fix_mm <- model.matrix(~ catchReg + month, catch)
# model matrix for predictions
mm_pred <- fix_mm %>% 
  unique() 
ord_mat <- mm_pred %>% 
  apply(., 1, sum) %>%
  sort() %>% 
  names()
mm_pred <- mm_pred[ord_mat, ]

data <- list(y1_i = catch$cpue,
             X1_ij = fix_mm,
             factor1k_i = yr_vec,
             nk1 = length(unique(yr_vec)),
             X1_pred_ij = mm_pred
)

# Fit simple model to initialize tmb
m1 <- lm(log_cpue ~ catchReg + month, data = catch)

parameters = list(
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
  log_phi = log(1.5),
  logit_p = boot::logit(0.8),
  z1_k = rep(0, length(unique(yr_vec))),
  log_sigma_zk1 = log(0.25)
)

# Factor key of unique combinations to generate predictions
fac_key <- catch %>% 
  select(catchReg, month) %>% 
  distinct() %>% 
  mutate(facs = paste(catchReg, month, sep = "_"),
         facs_n = as.numeric(as.factor(facs)))


# FIT --------------------------------------------------------------------------
# Compile
compile(here::here("R", "southCoast", "tmb", "tweedie_cpue_1re.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", "tweedie_cpue_1re")))

obj <- MakeADFun(data, parameters, random = c("z1_k"), 
                 DLL = "tweedie_cpue_1re")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr


# PREDICTIONS ------------------------------------------------------------------

## Plot predictions
N <- nrow(catch)
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
  geom_point(data = catch %>% filter(!cpue == 0), 
             aes(x = as.factor(month), y = log(cpue)), shape = 21, 
             alpha = 0.2) +
  geom_pointrange(data = pred_ci, aes(x = as.factor(month), y = log_pred_est,
                                      ymin = log_pred_low, 
                                      ymax = log_pred_up)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~catchReg)

ggplot() +
  geom_point(data = catch, aes(x = as.factor(month), y = cpue), shape = 21, 
             alpha = 0.2) +
  geom_pointrange(data = pred_ci, aes(x = as.factor(month), y = pred_est,
                                      ymin = pred_low, ymax = pred_up)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~catchReg)
