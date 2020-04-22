## Combined model fits - CWT
# April 21, 2020
# Fit combined multinomial/nb model to stock composition and abundance data
# using CWT recoveries from RMIS; compare to model predictions with GSI comp
# data from equivalent strata (basically equivalent to fit_combined)

library(tidyverse)
library(TMB)
library(ggplot2)

# Import Catch -----------------------------------------------------------------
# comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
#                            "reg3RollUpCatchProb.RDS")) 
comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                           "cwt_recovery_clean.rds")) %>% 
  rename(regName = pst_agg)

# month range dictated by ecological scale
month_range = c(1, 12)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                            "dailyCatch_WCVI.rds")) %>% 
  mutate(reg = factor(catchReg),
         area = as.factor(area),
         month_n = as.numeric(month),
         month = as.factor(month_n),
         year = as.factor(year)
  ) %>% 
  filter(!is.na(cpue),
         #first constrain by range
         !month_n < month_range[1],
         !month_n > month_range[2]
  ) %>% 
  droplevels() %>% 
  mutate(eff_z = as.numeric(scale(boatDays)),
         eff_z2 = eff_z^2
  ) %>% 
  arrange(catchReg, area, month)


# Clean Composition Data -------------------------------------------------------
source(here::here("R", "functions", "clean_composition_dat.R"))
comp_wide <- clean_comp(comp, month_range = month_range)

# Prep and Fit Model -----------------------------------------------------------

# make fixed effects factor key based on stock composition data
fac_dat <- comp_wide %>% 
  mutate(facs = as.factor(paste(as.character(catchReg), 
                                as.character(month), sep = "_")),
         facs = fct_relevel(facs, "NWVI_1", "NWVI_2", "NWVI_3", "NWVI_4", 
                            "NWVI_5", "NWVI_6",  "NWVI_7", "NWVI_8", "NWVI_9", 
                            "NWVI_10", "NWVI_11", "NWVI_12", "SWVI_1",  
                            "SWVI_2", "SWVI_3",  "SWVI_4", "SWVI_5", "SWVI_6",  
                            "SWVI_7", "SWVI_8",  "SWVI_9", "SWVI_10", 
                            "SWVI_11", "SWVI_12"),
         facs_n = as.numeric(facs) - 1) %>% 
  select(catchReg, month, facs, facs_n)
fac_key_eff <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n) %>%
  mutate(eff_z = mean(catch$eff_z),
         eff_z2 = mean(catch$eff_z2))

source(here::here("R", "functions", "prep_tmb_dat.R"))
inputs <- tmb_dat(catch, comp_wide, fac_dat, fac_key = fac_key_eff, mod = "nb")

## Make a function object
compile(here::here("R", "southCoast", "tmb", "nb_multinomial_1re.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", 
                           "nb_multinomial_1re")))
obj <- MakeADFun(inputs$data, inputs$parameters, random = c("z1_k", "z2_k"), 
                 DLL = "nb_multinomial_1re")

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr

# PREDICTIONS ------------------------------------------------------------------

log_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] #log pred of abundance
logit_probs <- ssdr[rownames(ssdr) %in% "logit_pred_prob", ] #logit probs of each category
pred_abund <- ssdr[rownames(ssdr) %in% "pred_abund_mg", ] #pred abundance of each category

comp_trim <- comp %>% 
  filter(!month_n < month_range[1],
         !month_n > month_range[2]) %>% 
  droplevels() %>% 
  dplyr::select(id, statArea, year, month, month_n, season, 
                regName, pres, catchReg)
n_groups <- length(unique(comp_trim$regName))

glimpse(inputs$data$X1_pred_ij)

pred_ci <- data.frame(stock = as.character(rep(unique(comp_trim$regName),
                                               each = 
                                                 length(unique(fac_key$facs_n)))), 
                      logit_prob_est = logit_probs[ , "Estimate"],
                      logit_prob_se =  logit_probs[ , "Std. Error"]) %>%
  mutate(facs_n = rep(fac_key$facs_n, times = n_groups)) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_low = plogis(logit_prob_est +
                                  (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est +
                                 (qnorm(0.975) * logit_prob_se)),
         abund_est = pred_abund[ , "Estimate"],
         abund_se =  pred_abund[ , "Std. Error"],
         abund_low = abund_est + (qnorm(0.025) * abund_se),
         abund_up = abund_est + (qnorm(0.975) * abund_se)) %>%
  left_join(., fac_key, by = c("facs_n")) 

# calculate raw summary data for comparison
raw_prop <- comp_trim %>% 
  group_by(catchReg, month, year, regName) %>%
  summarize(samp_g = length(unique(id))) %>% 
  group_by(catchReg, month, year) %>%
  mutate(samp_total = sum(samp_g)) %>% 
  ungroup() %>% 
  mutate(samp_g_ppn = samp_g / samp_total) %>% 
  rename(stock = regName) 

raw_abund <- catch %>% 
  group_by(catchReg, month, year) %>%
  summarize(sum_catch = sum(catch),
            sum_effort = sum(boatDays),
            agg_cpue = sum_catch / sum_effort) %>% 
  ungroup() %>% 
  left_join(., raw_prop, by = c("catchReg", "month", "year")) %>% 
  mutate(catch_g = samp_g_ppn * sum_catch,
         cpue_g = samp_g_ppn * agg_cpue,
         catchReg = as.factor(catchReg)) %>% 
  filter(!is.na(stock))

# focus on subset of stocks
stk_subset <- c("FR-early", "FR-late", "PSD", "CR-tule", "CR-bright", "WCVI")
plot_list <- map(list(pred_ci, raw_prop, raw_abund), function(x) 
  x %>% filter(stock %in% stk_subset)
  )

# combined estimates of stock-specific CPUE
# note that effort is standardized differently between raw data and predictions
ggplot() +
  geom_pointrange(data = plot_list[[1]], aes(x = month, y = abund_est,
                                      ymin = abund_low,
                                      ymax = abund_up), col = "red") +
  geom_point(data = plot_list[[3]], aes(x = month, y = cpue_g),
             alpha = 0.4) +
  facet_wrap(stock ~ catchReg, nrow = n_groups, scales = "free_y") +
  ggsidekick::theme_sleek()

# estimates of stock compostion
ggplot() +
  geom_pointrange(data = plot_list[[1]], aes(x = month, y = pred_prob, 
                                      ymin = pred_prob_low,
                                      ymax = pred_prob_up),
                  col = "red") +
  geom_point(data = plot_list[[2]],
             aes(x = month, y = samp_g_ppn),
             alpha = 0.4) +
  facet_wrap(stock ~ catchReg, nrow = n_groups, scales = "free_y") +
  ggsidekick::theme_sleek()
