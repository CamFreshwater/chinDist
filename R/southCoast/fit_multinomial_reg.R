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
# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                          "recIndProbsLong.rds"))
# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                          "wcviIndProbsLong.rds"))


# helper function to calculate aggregate probs
calc_max_prob <- function(grouped_data, full_data, thresh = 0.75) {
  out <- grouped_data %>% 
    summarize(agg_prob = sum(adj_prob)) %>% 
    arrange(id, desc(agg_prob)) %>%
    ungroup() %>% 
    group_by(id) %>% 
    mutate(max_assignment = max(agg_prob)) %>% 
    # Remove samples where top stock ID is less than threshold probability
    filter(!agg_prob < max_assignment, 
           !max_assignment < thresh) %>% 
    ungroup() %>% 
    distinct() %>% 
    left_join(.,
              full_data %>% 
                select(id:area_n) %>% 
                distinct(),
              by = "id") %>% 
    select(-max_assignment)
  
  colnames(out)[2] <- "agg"
  
  return(out)
}

# helper function to pool non-focal stocks
pool_aggs <- function(full_data) {
  full_data %>% 
    mutate(
      pst_agg = case_when(
        grepl("CR-", pst_agg) ~ "CR",
        grepl("CST", pst_agg) ~ "CA/OR/WA",
        pst_agg %in% c("NBC_SEAK", "WCVI") ~ "BC-coast",
        pst_agg %in% c("PSD", "SOG") ~ "SalSea",
        TRUE ~ pst_agg
      )
    )
}

rec_pst <- rec %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., rec)
comm_pst <- comm %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., comm)


# Trim based on dataset
trim_gen <- function(dat, month_range = c(1, 12)) {
  dat %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(id, region, area, year, month, season, agg, pres) 
}

gsi_trim <- trim_gen(comm_pst, month_range = c(4, 10))
gsi_trim <- trim_gen(rec_pst, month_range = c(6, 9))

table(gsi_trim$agg, gsi_trim$month)
table(gsi_trim$agg, gsi_trim$season, gsi_trim$region)
table(gsi_trim$agg, gsi_trim$month, gsi_trim$region)
table(gsi_trim$statArea, gsi_trim$month)
table(gsi_trim$month, gsi_trim$region)

# dummy dataset to replace missing values 
dum <- expand.grid(
  # year = unique(gsi_trim$year),
  month = unique(gsi_trim$month),
  # statArea = unique(gsi_trim$statArea),
  region = unique(gsi_trim$region),
  agg = unique(gsi_trim$agg),
  pres = 1)
#add random subset of yrs to avoid overparameterizing model
rand_yrs <- sample(unique(gsi_trim$year), size = nrow(dum), replace = TRUE)
dum$year <- rand_yrs
# add one value from each strata that was observed
# dum2 <- gsi_trim %>% 
#   select(month, region, agg, pres) %>% 
#   distinct() %>% 
#   mutate(year = sample(unique(gsi_trim$year), size = nrow(.), replace = TRUE)) %>% 
#   rbind(., dum)

gsi_wide <- gsi_trim %>%
  full_join(., dum, by = c("year", "month", "agg", "pres", "region")) %>%
  mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
  pivot_wider(., names_from = agg, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  droplevels() 

table(gsi_wide$agg, gsi_wide$month, gsi_wide$region)

# table(gsi_wide$year, gsi_wide$month, gsi_wide$agg, gsi_wide$statArea)
# table(gsi_wide$month, gsi_wide$agg, gsi_wide$statArea)


## RUN MODEL -------------------------------------------------------------------
compile(here::here("src", "multinomial_hier.cpp"))
dyn.load(dynlib(here::here("src", "multinomial_hier")))

# Prep data to pass to model
y_obs <- gsi_wide %>% 
  select(-c(id:dummy_id)) %>% 
  # select(Colmb, FrsrR, PgtSn, Other) %>% 
  as.matrix()
head(y_obs)

yr_vec <- as.numeric(gsi_wide$year) - 1
# fix_mm <- model.matrix(~ statArea + month, gsi_wide) #fixed covariates only
fix_mm <- model.matrix(~ region + month, gsi_wide) #fixed covariates only

#make combined factor levels (necessary for increasing speed of prob. estimates)
fac_dat <- gsi_wide %>% 
  mutate(facs = as.factor(paste(region, as.numeric(month), sep = "_")),
         #facs = as.factor(paste(statArea, month, year, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
  select(region, month, facs, facs_n)
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
  select(#-year,
    -facs, -facs_n) %>% 
  distinct()
# saveRDS(pred_ci_trim,
#         here::here("data", "gsiCatchData", "commTroll", 
#                    "multinomial_preds.RDS"))


# calculate raw proportion data for comparison
raw_prop <- gsi_trim %>% 
  group_by(region, month, year, agg) %>%
  summarize(samp_g = length(unique(id))) %>% 
  group_by(region, month, year) %>%
  mutate(samp_total = sum(samp_g)) %>% 
  ungroup() %>% 
  mutate(samp_g_ppn = samp_g / samp_total,
         stock = fct_reorder(agg, desc(samp_g_ppn))) 


ggplot(pred_ci_trim) +
  geom_point(data = raw_prop, 
             aes(x = month, y = samp_g_ppn, fill = region),
             shape = 21, alpha = 0.4, position = position_dodge(0.6)) +
  geom_pointrange(aes(x = as.factor(month), y = pred_prob, ymin = pred_prob_low,
                      ymax = pred_prob_up, fill = region), 
                  shape = 21, size = 0.4) +
  labs(y = "Probability", x = "Month") +
  facet_wrap(~ stock) +
  ggsidekick::theme_sleek()





ggplot() +
  geom_point(data = raw_prop, 
             aes(x = month, y = samp_g_ppn, fill = catchReg),
             shape = 21, alpha = 0.3, position = position_dodge(0.6)) +
  geom_pointrange(data = pred_ci, 
                  aes(x = month, y = pred_prob, ymin = pred_prob_low, 
                      ymax = pred_prob_up, fill = catchReg), 
                  shape = 21, position = position_dodge(0.6)) +
  facet_wrap(~stock, ncol = 2, scales = "free_y") +
  scale_fill_viridis_d(option = "C", name = "Catch Region") +
  labs(x = "Month", y = "Predicted Encounter Probability") +
  ggsidekick::theme_sleek() +
  theme(legend.position="top")