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
  filter(!statArea %in% c("24", "26")) %>% 
  droplevels()
table(reg3_trim$regName, reg3_trim$month)
table(reg3_trim$regName, reg3_trim$season, reg3_trim$catchReg)
table(reg3_trim$statArea, reg3_trim$month, reg3_trim$regName)
table(reg3_trim$statArea, reg3_trim$catchReg)

# dummy dataset to replace missing values 
dum <- expand.grid(
  # year = unique(reg3_trim$year),
  month = unique(reg3_trim$month),
  statArea = unique(reg3_trim$statArea),
  # catchReg = unique(reg3_trim$catchReg),
  regName = unique(reg3_trim$regName),
  pres = 1)
#add random subset of yrs to avoid overparameterizing model
rand_yrs <- sample(unique(reg3_trim$year), size = nrow(dum), replace = TRUE)
dum$year <- rand_yrs
# add one value from each strata that was observed
dum2 <- reg3_trim %>% 
  select(month, statArea, regName, pres) %>% 
  distinct() %>% 
  mutate(year = sample(unique(reg3_trim$year), size = nrow(.), replace = TRUE)) %>% 
  rbind(., dum)

gsi_wide <- reg3_trim %>%
  # full_join(., dum, by = c("year", "month", "regName", "pres", "catchReg")) %>%
  full_join(., dum2, by = c("year", "month", "regName", "pres", "statArea")) %>%
  mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  droplevels()

# alternative way of avoiding zero values, add 1 to each category
# dum <- expand.grid(
#   statArea = unique(reg3_trim$statArea),
#   year = unique(reg3_trim$year),
#   month = unique(reg3_trim$month),
#   regName = unique(reg3_trim$regName),
#   pres = 1)
# 
# gsi_wide <- reg3_trim %>% 
#   select(-flatFileID, -catchReg, -season) %>% 
#   rbind(., dum) %>% 
#   mutate(id = seq(from = 1, to = nrow(.), by = 1)) %>% 
#   pivot_wider(., names_from = regName, values_from = pres) %>%
#   mutate_if(is.numeric, ~replace_na(., 0)) %>% 
#   droplevels() %>% 
#   arrange(statArea, year, month)

table(gsi_wide$year, gsi_wide$month, gsi_wide$regName, gsi_wide$statArea)
table(gsi_wide$month, gsi_wide$regName, gsi_wide$statArea)


## RUN MODEL -------------------------------------------------------------------
compile(here::here("R", "southCoast", "tmb", "multinomial_hier.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", "multinomial_hier")))

# Prep data to pass to model
y_obs <- gsi_wide %>% 
  select(Colmb, FrsrR, PgtSn, Other) %>% 
  as.matrix()
yr_vec <- as.numeric(gsi_wide$year) - 1
fix_mm <- model.matrix(~ statArea + month, gsi_wide) #fixed covariates only

#make combined factor levels (necessary for increasing speed of prob. estimates)
fac_dat <- gsi_wide %>% 
  mutate(facs = as.factor(paste(statArea, month, year, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
  select(statArea, month, year, facs, facs_n)
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
sdr2 <- sdreport(obj, getReportCovariance = TRUE)
sdr
sdr2

ssdr <- summary(sdr)
ssdr

## Plot predictions
k <- ncol(y_obs) # number of stocks
stk_names <- colnames(y_obs)
N <- nrow(y_obs)

#fixed effects
logit_probs_mat <- ssdr[rownames(ssdr) %in% "logit_probs_out", ]
stock_n_vec <- as.character(rep(1:k, each = length(unique(fac_key$facs_n))))
stock_vec <- as.character(rep(stk_names, each = length(unique(fac_key$facs_n))))
pred_ci_pool <- data.frame(stock_n = stock_n_vec, 
                           stock = stock_vec, 
                           logit_prob_est = logit_probs_mat[ , "Estimate"],
                           logit_prob_se =  logit_probs_mat[ , "Std. Error"]) %>%
  mutate(ests = "pool",
         facs_n = rep(fac_key$facs_n, times = k))

logit_probs_mat_fe <- ssdr[rownames(ssdr) %in% "logit_probs_out_fe", ]
pred_ci <- data.frame(stock_n = stock_n_vec, 
                      stock = stock_vec, 
                      logit_prob_est = logit_probs_mat_fe[ , "Estimate"],
                      logit_prob_se =  logit_probs_mat_fe[ , "Std. Error"]) %>%
  mutate(ests = "fix",
         facs_n = rep(fac_key$facs_n, times = k)) %>% 
  rbind(pred_ci_pool, .) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_se = plogis(logit_prob_se),
         pred_prob_low = plogis(logit_prob_est +
                                  (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est +
                                 (qnorm(0.975) * logit_prob_se))) %>%
  left_join(., fac_key, by = "facs_n") %>%
  select(-logit_prob_est, -logit_prob_se) 

pred_ci_trim <- pred_ci %>% 
  filter(ests == "fix") %>% 
  select(-year, -facs, -facs_n) %>% 
  distinct()
saveRDS(pred_ci_trim,
        here::here("data", "gsiCatchData", "commTroll", 
                   "multinomial_preds.RDS"))

ggplot(pred_ci_trim) +
  geom_pointrange(aes(x = as.factor(month), y = pred_prob, ymin = pred_prob_low, 
                      ymax = pred_prob_up, fill = statArea), shape = 21) +
  labs(y = "Probability", x = "Month") +
  facet_wrap(~ stock)

