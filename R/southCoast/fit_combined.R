## Combined model fits
# April 3, 2020
# Fit combined multionomial/tweedie or multinomial/nb model to 
# stock composition and abundance data
# aggregate probability reflects summed probabilities of a given region of 
# origin (ie reg1 or 3) for a given individual

library(tidyverse)
library(TMB)
library(ggplot2)


# Import Catch -----------------------------------------------------------------
gsi <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                          "reg1RollUpCatchProb_FraserB.RDS"))
# gsi <- readRDS(here::here("data", "gsiCatchData", "commTroll",
#                           "reg3RollUpCatchProb.RDS"))

# pull months that are common to both strata and subset
comm_months <- gsi %>% 
  select(catchReg, month) %>% 
  distinct() %>% 
  split(., .$catchReg) %>% 
  map(., function(x) x %>% pull(as.numeric(month))) %>% 
  Reduce(intersect, .)
#select a range of months to utilize based on fraser focus or aggregate 
#and spatial strata based on data availability
if(is.null(gsi$pscName)) {
  month_range <- c(1, 12)
  sp_var = "area"
} else {
  month_range <- c(3, 10)
  sp_var = "reg"
}
                      
#add dummy catch data for one month that's missing based on gsi data
min_catch <- gsi %>% 
  filter(catchReg == "SWVI",
         month == "7") %>%
  group_by(year) %>% 
  summarise(n = length(unique(flatFileID)),
            n_days = length(unique(jDay)) * 3)
min_catch_dat <- data.frame(catch = min_catch$n,
                            catchReg = "SWVI",
                            area = NA,
                            month = "7",
                            jDay = NA,
                            year = min_catch$year,
                            boatDays = min_catch$n_days) %>% 
  mutate(cpue = catch / boatDays) %>% 
  select(catchReg:year, catch, boatDays, cpue)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                            "dailyCatch_WCVI.rds")) %>% 
  rbind(., min_catch_dat) %>% 
  mutate(reg = factor(catchReg),
         area = as.factor(area),
         month_n = as.numeric(month),
         month = as.factor(month_n),
         year = as.factor(year)
  ) %>% 
  filter(!is.na(cpue),
         #first constrain by range
         !month_n < month_range[1],
         !month_n > month_range[2],
         #then drop missing months
         month_n %in% comm_months) %>% 
  droplevels() %>% 
  mutate(eff_z = as.numeric(scale(boatDays)),
         eff_z2 = eff_z^2
  ) %>% 
  arrange(catchReg, area, month)

# if(is.null(gsi$pscName)) {
#   catch <- catch %>% 
#     mutate(sp_var = "area")
# } else {
#   month_range <- c(4, 10)
#   sp_var = "reg"
# }


# Import Genetics --------------------------------------------------------------
# Trim based on dataset
trim_gen <- function(dat, month_range = c(1, 12)) {
  dat %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(flatFileID, statArea, year, month, month_n, season, 
                  regName, pres, catchReg) 
}

gsi_trim <- trim_gen(gsi, month_range = month_range)

# sp_var_g <- gsi_trim %>% pull(sp_var)
table(gsi_trim$regName, gsi_trim$month, gsi_trim$catchReg)
# sp_var_c <- catch %>% pull(sp_var)
table(catch$month, catch$catchReg)

# multinomial model predictions are a little wonky w/ missing values - this is 
# one infilling option
# dummy dataset to replace missing values 
stock_names <- unique(gsi_trim$regName)
#retain only values from the dummy set that are not already present
missing_sample_events <- expand.grid(
  flatFileID = NA,
  statArea = NA,
  year = NA,
  month = unique(gsi_trim$month),
  month_n = NA,
  season = NA,
  regName = stock_names,
  pres = 1,
  catchReg = NA#unique(gsi_trim$catchReg)
  ) %>%
  anti_join(., gsi_trim,
                      by = c("month", "regName", "pres"#, "catchReg"
                             )) %>%
  select(-regName) %>%
  distinct()
#add random subset of yrs to avoid overparameterizing model
rand_yrs <- sample(unique(gsi_trim$year), size = nrow(missing_sample_events),
                   replace = TRUE)
rand_catch_reg <- sample(unique(gsi_trim$catchReg), 
                         size = nrow(missing_sample_events),
                   replace = TRUE)
missing_sample_events$year <- rand_yrs
missing_sample_events$catchReg <- rand_catch_reg

#duplicate for all stocks (to balance the addition)
infill_data <- do.call("rbind", replicate(length(stock_names),
                                          missing_sample_events,
                                          simplify = FALSE)) %>%
  mutate(regName = rep(stock_names, each = nrow(missing_sample_events))) %>%
  select(flatFileID:season, regName, pres, catchReg)

# combine and check for no zeros
temp <- gsi_trim %>%
  rbind(., infill_data)
table(temp$regName, temp$month, temp$catchReg)
table(catch$month, catch$catchReg)

# spread
gsi_wide <- temp %>% 
  mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  droplevels() 


# Prep Data for TMB ------------------------------------------------------------

# fixed effects model matrices that either include or exclude effort
# form_c <- paste("~ -1 + (catchReg : month)")
sp_t_ch <- paste("~ catchReg + month")
sp_t_form <- formula(sp_t_ch)
sp_t_eff_ch <- paste("~ catchReg + month + eff_z + eff_z2")
sp_t_eff_form <- formula(sp_t_eff_ch)

fix_mm_gsi <- model.matrix(sp_t_form, data = gsi_wide)
fix_mm_c <- if (mod == "tweedie") {
  model.matrix(sp_t_form, data = catch)
} else {
  model.matrix(sp_t_eff_form, data = catch)
} 

#helper function to convert factors 
fct_to_tmb_num <- function(x) {
  as.numeric(as.factor(as.character(x))) - 1
}

# fit dummy model to speed up tweedie estimatess
form_dum <- if (mod == "tweedie") {
  formula(paste("log(cpue + 0.0001)", sp_t_ch, sep = ""))
} else {
  formula(paste("log(catch + 0.0001)", sp_t_eff_ch, sep = ""))
}
m1 <- lm(form_dum, data = catch)

# make fixed effects factor key based on stock composition data
fac_dat <- gsi_wide %>% 
  mutate(facs = as.factor(paste(as.character(catchReg), 
                                as.character(month), sep = "_")),
         facs = fct_reorder2(facs, 
                             desc(as.numeric(as.character(month))),
                             desc(catchReg)),
         facs = fct_relevel(facs, "SWVI_10", after = 8),
         facs = fct_relevel(facs, "NWVI_10", after = Inf),
         facs_n = as.numeric(facs) - 1) %>% 
  select(catchReg, month, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)

# predictive model matrix for abundance based on fac key with the potential for
# the addition of covariate effects
fac_key_eff <- fac_key %>% 
  mutate(eff_z = mean(catch$eff_z),
         eff_z2 = mean(catch$eff_z2))
mm_pred <- if (mod == "tweedie") {
  model.matrix(sp_t_form, data = fac_key)
} else {
  model.matrix(sp_t_eff_form, data = fac_key_eff)
}

if(!identical(colnames(mm_pred), colnames(fix_mm_c))) {
  warning("Mismatch between estimated and predicted abundance model matrices.")
}

# observed stock composition
y_obs <- gsi_wide %>% 
  select(-c(flatFileID:dummy_id)) %>% 
  as.matrix()
head(y_obs)

# other parameters
n_groups <- ncol(y_obs)
b1_n <- length(coef(m1))
b2_n <- ncol(mm_pred) * (n_groups - 1) #each par est for each non-ref level
# vectors of random effects
fac1k <- fct_to_tmb_num(catch$year)
fac2k <- fct_to_tmb_num(gsi_wide$year)
nk1 <- length(unique(fac1k))
nk2 <- length(unique(fac2k))

response <- if (mod == "tweedie") {
  catch$cpue
} else{
  catch$catch
}

# combine
data <- list(
  #abundance data
  y1_i = response,
  X1_ij = fix_mm_c,
  factor1k_i = fac1k,
  nk1 = nk1,
  X1_pred_ij = mm_pred,
  #composition data
  y2_ig = y_obs,
  X2_ij = fix_mm_gsi,
  factor2k_i = fac2k,
  nk2 = nk2,
  m2_all_fac = fac_dat$facs_n,
  m2_fac_key = fac_key$facs_n
)

parameters = list(
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
  log_phi = log(1.1),
  logit_p = boot::logit(0.8),
  z1_k = rep(0, length(unique(fac1k))),
  log_sigma_zk1 = log(0.25),
  b2_jg = matrix(0, nrow = ncol(fix_mm_gsi), ncol = (n_groups - 1)), 
  z2_k = rep(0, times = length(unique(fac2k))),
  log_sigma_zk2 = log(0.25)
)

if(!mod == "tweedie") {
  parameters[["logit_p"]] <- NULL
}

# Fit Model --------------------------------------------------------------------

## Make a function object
compile(here::here("R", "southCoast", "tmb", "tweedie_multinomial_1re.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", 
                           "tweedie_multinomial_1re")))
obj <- MakeADFun(data, parameters, random = c("z1_k", "z2_k"), 
                 DLL = "tweedie_multinomial_1re")

compile(here::here("R", "southCoast", "tmb", "nb_multinomial_1re.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", 
                           "nb_multinomial_1re")))
obj <- MakeADFun(data, parameters, random = c("z1_k", "z2_k"), 
                 DLL = "nb_multinomial_1re")

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr

# saveRDS(ssdr, here::here("generatedData", "model_fits", "twmult_ssdr_agg.RDS"))
# saveRDS(ssdr, here::here("generatedData", "model_fits", "twmult_ssdr_frB.RDS"))

# PREDICTIONS ------------------------------------------------------------------

#frB versions have a different subset of stocks and months at finer spatial scale
# ssdr <- readRDS(here::here("generatedData", "model_fits", "twmult_ssdr_frB.RDS"))
# ssdr <- readRDS(here::here("generatedData", "model_fits", "twmult_ssdr_agg.RDS"))


log_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] #log pred of abundance
logit_probs <- ssdr[rownames(ssdr) %in% "logit_pred_prob", ] #logit probs of each category
pred_abund <- ssdr[rownames(ssdr) %in% "pred_abund_mg", ] #pred abundance of each category

pred_ci <- data.frame(stock = as.character(rep(unique(gsi_trim$regName),
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
raw_prop <- gsi_trim %>% 
  group_by(catchReg, month, year, regName) %>%
  summarize(samp_g = length(unique(flatFileID))) %>% 
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
  mutate(catch_g = samp_g_ppn * z_catch,
         cpue_g = samp_g_ppn * agg_cpue,
         catchReg = as.factor(catchReg)) %>% 
  # rename(stock = regName) %>% 
  filter(!is.na(stock))

agg_abund <- data.frame(
  facs_n = fac_key$facs_n,
  raw_abund_est = log_pred[ , "Estimate"],
  raw_abund_se = log_pred[ , "Std. Error"]) %>% 
  mutate(
    raw_mu = exp(raw_abund_est),
    raw_abund_low = exp(raw_abund_est + (qnorm(0.025) * raw_abund_se)),
    raw_abund_up = exp(raw_abund_est + (qnorm(0.975) * raw_abund_se))
  ) %>% 
  left_join(pred_ci, ., by = "facs_n")


# combined estimates of stock-specific CPUE
ggplot() +
  geom_pointrange(data = pred_ci, aes(x = month, y = abund_est,
                                      ymin = abund_low,
                                      ymax = abund_up), col = "red") +
  # geom_point(data = raw_abund, aes(x = month, y = catch_g),
  #            alpha = 0.4) +
  geom_point(data = raw_abund, aes(x = month, y = cpue_g),
             alpha = 0.4) +
  facet_wrap(stock ~ catchReg, nrow = n_groups, scales = "free_y") +
  ggsidekick::theme_sleek()

# estimates of stock compostion
ggplot() +
  geom_pointrange(data = pred_ci, aes(x = month, y = pred_prob, 
                                      ymin = pred_prob_low,
                                      ymax = pred_prob_up),
                  col = "red") +
  geom_point(data = raw_prop,
             aes(x = month, y = samp_g_ppn),
             alpha = 0.4) +
  facet_wrap(stock ~ catchReg, nrow = n_groups, scales = "free_y") +
  ggsidekick::theme_sleek()

# estimates of aggregate CPUE
ggplot() +
  geom_boxplot(data = catch, aes(x = month, y = cpue),  
             alpha = 0.4) +
  geom_pointrange(data = agg_abund, aes(x =  month, y = raw_mu,
                                  ymin = raw_abund_low,
                                  ymax = raw_abund_up), color= "red") +
  facet_wrap(~ catchReg, nrow = 2, scales = "free_y") +
  ggsidekick::theme_sleek()

## export plotting data for Rmd
list(catch = catch, pred_ci = pred_ci, raw_prop = raw_prop, raw_abund = raw_abund, 
           raw_agg_abund = agg_abund) %>% 
  saveRDS(., here::here("generatedData", "model_fits", "frB_plot_list.RDS"))


## estimates of effort effects on catch
n_betas <- length(coef(m1))
abund_b <- ssdr[rownames(ssdr) %in% "b1_j", ]
eff_b <- abund_b[c(1, n_betas - 1, n_betas) , 1]

ggplot(catch) +
  geom_point(aes(x = eff_z, y = catch), 
             alpha = 0.2) +
  # lims(x = c(0, 4)) +
  stat_function(fun = function(x) exp(eff_b[1] + eff_b[2]*x + eff_b[3]*x^2)) 


# look at predicted cumulative abundance
plot_cum_dens <- function(dat, station = "JDF") {
  dat %>% 
    filter(project_name == station) %>% 
    ggplot(., aes(yday, colour = CU)) +
    stat_ecdf()
}

pred_ci %>% 
  arrange(month) %>% 
  group_by(stock, catchReg) %>% 
  mutate(total_abund = sum(abund_est),
         cum_abund = cumsum(abund_est) / total_abund) %>%
  ungroup() %>% 
  ggplot(., aes(x = month, y = cum_abund, colour = stock)) +
  geom_line() +
  geom_point() +
  facet_wrap(~catchReg)

