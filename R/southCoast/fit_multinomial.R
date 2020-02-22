## Multinomial model fit
# January 20, 2020
# Fit multinomial model to stock composition data 
# (assumes perfect GSI assignment)

library(tidyverse)
library(TMB)
library(ggplot2)


reg3_long <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "reg3RollUpCatchProb.RDS")) %>% 
  #remove all secondary IDs and uncertain IDs (less than 80%)
  group_by(flatFileID) %>% 
  mutate(maxProb = max(aggProb)) %>% 
  filter(!aggProb != maxProb, 
         !maxProb < 0.8) %>% 
  ungroup() %>% 
  # aggregate as necessary
  mutate(
    month = as.factor(month),
    year =  as.factor(year),
    regName = case_when(
      regName == "Oregon/California" ~ "Or_Cl",
      regName %in% c("Columbia", "Snake") ~ "Columbia",
      regName %in% c("Coastal Washington", "Washington Coast", 
                     "Alaska South SE", "North/Central BC", "SOG") ~ "Other",
      TRUE ~ regName
    ),
    regName = as.factor(abbreviate(regName, minlength = 5)),
    pres = 1,
    area_n = as.numeric(as.character(statArea)),
    catchReg = case_when(
      area_n < 125 & area_n > 27 ~ "SWVI",
      area_n < 25 ~ "SWVI",
      TRUE ~ "NWVI"
    )) %>%
  dplyr::select(flatFileID, statArea, year, month, regName, pres, catchReg)

reg3_trim <- reg3_long %>% 
  filter(statArea == "123") %>% 
  droplevels()
table(reg3_trim$regName, reg3_trim$month)

# dummy dataset to replace missing values 
dum <- expand.grid(#year = unique(reg3_long$year), 
  month = unique(reg3_trim$month),
  regName = unique(reg3_trim$regName), 
  pres = 1)

gsi_wide <- reg3_trim %>% 
  full_join(., dum, by = c("month", "regName", "pres")) %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0))

# Prep data to pass to model
obs_mat <- gsi_wide %>% 
  select(FrsrR:WCVI) %>% 
  as.matrix()


.X <- model.matrix(~ month, gsi_wide)

data <- list(cov = .X, y_obs = obs_mat)


## RUN MODEL -------------------------------------------------------------------
# compile("C:/github/juvenile-salmon-index/R/multinomialPractice/multinomial_generic.cpp")
# dyn.load(dynlib("C:/github/juvenile-salmon-index/R/multinomialPractice/multinomial_generic"))

compile(here::here("R", "southCoast", "multinomial_generic.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "multinomial_generic")))

pars <- list(betas = matrix(data = 0, nrow = ncol(.X), 
                            ncol = ncol(obs_mat) - 1))

## Make a function object
obj <- MakeADFun(data, pars, DLL="multinomial_generic")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

ssdr <- summary(sdr)
ssdr

r <- obj$report()
r$probs
r$log_odds
r$logit_probs


## Plot predictions
k <- ncol(obs_mat) # number of stocks
stk_names <- colnames(obs_mat)
N <- nrow(obs_mat)
months <- gsi_wide$month
years <- gsi_wide$year

logit_probs_mat <- ssdr[rownames(ssdr) %in% "logit_probs", ]
# logit_probs_mat <- r$logit_probs
pred_ci <- data.frame(stock = rep(stk_names, each = N),
                      month = rep(months, times = k),
                      #year = as.factor(rep(years, times = k)),
                      logit_prob_est = logit_probs_mat[ , "Estimate"],
                      logit_prob_se =  logit_probs_mat[ , "Std. Error"]) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_low = plogis(logit_prob_est + (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est + (qnorm(0.975) * logit_prob_se))
  ) %>%
  distinct()

ggplot(pred_ci) +
  # geom_point(aes(x = year, y = pred_prob)) +
  geom_pointrange(aes(x = month, y = pred_prob, ymin = pred_prob_low, 
                      ymax = pred_prob_up)) +
  # geom_ribbon(aes(x = jday, ymin = pred_prob_low, ymax = pred_prob_up), 
  #             fill = "#bfd3e6") +
  # geom_line(aes(x = jday, y = pred_prob), col = "#810f7c", size = 1) +
  facet_grid( ~ stock) +
  labs(y = "Probability", x = "Year")




