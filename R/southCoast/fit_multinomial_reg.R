## Multinomial model fit
# January 20, 2020
# Fit multinomial model to stock composition data at regional aggregate level
# (assumes perfect GSI assignment)

library(tidyverse)
library(TMB)
library(ggplot2)

reg3 <- readRDS(here::here("generatedData", "daily_gsi_maxprob.rds")) %>% 
  # aggregate as necessary
  mutate(
    season_c = case_when(
      month %in% c("12", "1", "2") ~ "w",
      month %in% c("3", "4", "5") ~ "sp",
      month %in% c("6", "7", "8") ~ "su",
      month %in% c("9", "10", "11") ~ "f"
    ),
    season = fct_relevel(season_c, "sp", "su", "f", "w"),
    month = as.factor(month),
    year =  as.factor(year),
    regName = as.character(regName),
    regName = case_when(
      # regName == "Oregon/California" ~ "Or_Cl",
      regName %in% c("Columbia", "Snake") ~ "Columbia",
      regName %in% c("Coastal Washington", "Washington Coast",
                     "Alaska South SE", "North/Central BC", "SOG", 
                     "Oregon/California", "ECVI", "WCVI") ~ "Other",
      # regName %in% c("Coastal Washington", "Washington Coast",
      #                "Alaska South SE", "North/Central BC", "SOG") ~ "Other",
      TRUE ~ regName
    ),
    regName = as.factor(abbreviate(regName, minlength = 5)),
    pres = 1,
    area_n = as.numeric(as.character(statArea)),
    catchReg = case_when(
      area_n < 125 & area_n > 27 ~ "SWVI",
      area_n < 25 ~ "SWVI",
      TRUE ~ "NWVI"
    ),
    catchReg = as.factor(catchReg), 
    statArea = as.factor(statArea)) %>%
  dplyr::select(flatFileID, statArea, year, month, season, regName, pres, 
                catchReg)

reg3_trim <- reg3 %>% 
  droplevels()
table(reg3_trim$regName, reg3_trim$month)
table(reg3_trim$regName, reg3_trim$season, reg3_trim$catchReg)
table(reg3_trim$regName, reg3_trim$year)

# dummy dataset to replace missing values 
dum <- expand.grid(
  year = unique(reg3_trim$year),
  # season = unique(reg3_trim$season),
  month = unique(reg3_trim$month),
  catchReg = unique(reg3_trim$catchReg),
  regName = unique(reg3_trim$regName),
  pres = 1)

gsi_wide <- reg3_trim %>% 
  # sample_n(2000) %>%
  full_join(., dum, by = c("year", "month", "regName", "pres", "catchReg")) %>%
  # arrange(regName) %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  # filter(season %in% c("sp", "su")) %>% 
  droplevels() 

table(gsi_wide$year, gsi_wide$month, gsi_wide$regName, gsi_wide$catchReg)


## RUN MODEL -------------------------------------------------------------------
compile(here::here("R", "southCoast", "tmb", "multinomial_hier.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", "multinomial_hier")))

# Prep data to pass to model
y_obs <- gsi_wide %>% 
  select(Colmb, FrsrR, PgtSn, Other) %>% 
  as.matrix()
yr_vec <- as.numeric(gsi_wide$year) - 1
fix_mm <- model.matrix(~ catchReg + month, gsi_wide) #fixed covariates only

#make combined factor levels (necessary for increasing speed of prob. estimates)
fac_dat <- gsi_wide %>% 
  mutate(facs = as.factor(paste(catchReg, month, year, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
  select(catchReg, month, year, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)

data <- list(y_obs = y_obs, #obs
             rfac = yr_vec, #random intercepts
             fx_cov = fix_mm, #fixed cov model matrix
             n_rfac = length(unique(yr_vec)), #number of random intercepts
             all_fac = fac_dat$facs_n, # vector of factor combinations
             fac_key = fac_key$facs_n #ordered unique factor combos in fac_vec
) 
parameters <- list(z_rfac = rep(0, times = length(unique(yr_vec))),
                   z_ints = matrix(0, nrow = ncol(fix_mm), 
                                   ncol = ncol(y_obs) - 1),
                   log_sigma_rfac = 0)

## Make a function object
obj <- MakeADFun(data, parameters, random = c("z_rfac"), 
                 DLL = "multinomial_hier")

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
k <- ncol(y_obs) # number of stocks
stk_names <- colnames(y_obs)
N <- nrow(y_obs)

#fixed effects
catch_reg <- gsi_wide$catchReg
years <- gsi_wide$year


logit_probs_mat <- ssdr[rownames(ssdr) %in% "logit_probs", ]
pred_ci <- data.frame(stock = rep(stk_names, each = N),
                      # month = rep(months, times = k),
                      catch_reg = rep(catch_reg, times = k),
                      # year = rep(years, times = k),
                      logit_prob_est = logit_probs_mat[ , "Estimate"],
                      logit_prob_se =  logit_probs_mat[ , "Std. Error"]) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_low = plogis(logit_prob_est + (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est + (qnorm(0.975) * logit_prob_se))
  ) %>%
  distinct()

ggplot(pred_ci) +
  # geom_point(aes(x = year, y = pred_prob)) +
  geom_pointrange(aes(x = stock, y = pred_prob, ymin = pred_prob_low, 
                      ymax = pred_prob_up)) +
  # geom_ribbon(aes(x = jday, ymin = pred_prob_low, ymax = pred_prob_up), 
  #             fill = "#bfd3e6") +
  # geom_line(aes(x = jday, y = pred_prob), col = "#810f7c", size = 1) +
  labs(y = "Probability", x = "Stock") +
  facet_wrap(~catch_reg)




