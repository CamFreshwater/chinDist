## Multinomial model fit
# January 20, 2020
# Fit multinomial model to stock composition data at various aggregate levels
# (assumes perfect GSI assignment); in both cases aggregate probability reflects
# summed probabilities of a given region of origin (ie reg1 or 3) for a given
# individual

library(tidyverse)
library(TMB)
library(ggplot2)


# Import various regional roll ups and initial cleaning
# region3 roll up
reg3 <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "reg3RollUpCatchProb.RDS")) 
# region1 roll up
reg1_fr <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "reg1RollUpCatchProb_Fraser.RDS")) 

# Trim based on dataset
trim_gen <- function(dat, month_range = c(1, 12)) {
  dat %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(flatFileID, statArea, year, month, season, regName, 
                  pres, catchReg) 
}

gsi_trim <- trim_gen(reg1_fr, month_range = c(4, 10))

table(gsi_trim$regName, gsi_trim$month)
table(gsi_trim$regName, gsi_trim$season, gsi_trim$catchReg)
table(gsi_trim$regName, gsi_trim$month, gsi_trim$catchReg)
table(gsi_trim$statArea, gsi_trim$month)
table(gsi_trim$month, gsi_trim$catchReg)

# dummy dataset to replace missing values 
dum <- expand.grid(
  # year = unique(gsi_trim$year),
  month = unique(gsi_trim$month),
  # statArea = unique(gsi_trim$statArea),
  catchReg = unique(gsi_trim$catchReg),
  regName = unique(gsi_trim$regName),
  pres = 1)
#add random subset of yrs to avoid overparameterizing model
rand_yrs <- sample(unique(gsi_trim$year), size = nrow(dum), replace = TRUE)
dum$year <- rand_yrs
# add one value from each strata that was observed
# dum2 <- gsi_trim %>% 
#   select(month, catchReg, regName, pres) %>% 
#   distinct() %>% 
#   mutate(year = sample(unique(gsi_trim$year), size = nrow(.), replace = TRUE)) %>% 
#   rbind(., dum)

gsi_wide <- gsi_trim %>%
  full_join(., dum, by = c("year", "month", "regName", "pres", "catchReg")) %>%
  mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  droplevels()

table(gsi_wide$regName, gsi_wide$month, gsi_wide$catchReg)

# table(gsi_wide$year, gsi_wide$month, gsi_wide$regName, gsi_wide$statArea)
# table(gsi_wide$month, gsi_wide$regName, gsi_wide$statArea)


## RUN MODEL -------------------------------------------------------------------
compile(here::here("R", "southCoast", "tmb", "multinomial_hier.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", "multinomial_hier")))

# Prep data to pass to model
y_obs <- gsi_wide %>% 
  select(-c(flatFileID:dummy_id)) %>% 
  # select(Colmb, FrsrR, PgtSn, Other) %>% 
  as.matrix()
head(y_obs)

yr_vec <- as.numeric(gsi_wide$year) - 1
# fix_mm <- model.matrix(~ statArea + month, gsi_wide) #fixed covariates only
fix_mm <- model.matrix(~ catchReg + month, gsi_wide) #fixed covariates only

#make combined factor levels (necessary for increasing speed of prob. estimates)
fac_dat <- gsi_wide %>% 
  mutate(facs = as.factor(paste(catchReg, month, year, sep = "_")),
         #facs = as.factor(paste(statArea, month, year, sep = "_")),
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
# saveRDS(pred_ci_trim,
#         here::here("data", "gsiCatchData", "commTroll", 
#                    "multinomial_preds.RDS"))

ggplot(pred_ci_trim) +
  geom_pointrange(aes(x = as.factor(month), y = pred_prob, ymin = pred_prob_low, 
                      ymax = pred_prob_up, fill = catchReg), shape = 21) +
  labs(y = "Probability", x = "Month") +
  facet_wrap(~ stock)

