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
# comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
#                           "reg1RollUpCatchProb_FraserB.RDS"))
comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                          "reg3RollUpCatchProb.RDS"))

# month range dictated by ecological scale
month_range = c(1, 10)

# pull months that are common to both strata and subset
# comm_months <- comp %>% 
#   select(catchReg, month) %>% 
#   distinct() %>% 
#   split(., .$catchReg) %>% 
#   map(., function(x) x %>% pull(as.numeric(month))) %>% 
#   Reduce(intersect, .)

#add dummy catch data for one month that's missing logbook data based on observed
# catch in comp data
min_catch <- comp %>% 
  filter(catchReg == "SWVI",
         month == "7") %>%
  group_by(year) %>% 
  summarise(n = length(unique(id)),
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
         !month_n > month_range[2]
         # ,
         # #then drop missing months
         # month_n %in% comm_months
         ) %>% 
  droplevels() %>% 
  mutate(eff_z = as.numeric(scale(boatDays)),
         eff_z2 = eff_z^2
  ) %>% 
  arrange(catchReg, area, month)


# Import Genetics --------------------------------------------------------------
# Function to trim and infill comp dataset
clean_gsi <- function(comp, month_range, check_tables = FALSE) {
  comp_trim <- comp %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(id, statArea, year, month, month_n, season, 
                  regName, pres, catchReg) 
  
  # multinomial model predictions are a little wonky w/ missing values - this is 
  # one infilling option
  # dummy dataset to replace missing values 
  stock_names <- unique(comp_trim$regName)
  #retain only values from the dummy set that are not already present
  missing_sample_events <- expand.grid(
    id = NA,
    statArea = NA,
    year = NA,
    month = unique(comp_trim$month),
    month_n = NA,
    season = NA,
    regName = stock_names,
    pres = 1,
    catchReg = NA
  ) %>%
    anti_join(., comp_trim,
              by = c("month", "regName", "pres")
              ) %>%
    select(-regName) %>%
    distinct()
  #add random subset of yrs to avoid overparameterizing model
  rand_yrs <- sample(unique(comp_trim$year), size = nrow(missing_sample_events),
                     replace = TRUE)
  rand_catch_reg <- sample(unique(comp_trim$catchReg), 
                           size = nrow(missing_sample_events),
                           replace = TRUE)
  missing_sample_events$year <- rand_yrs
  missing_sample_events$catchReg <- rand_catch_reg
  
  #duplicate for all stocks (to balance the addition)
  infill_data <- do.call("rbind", replicate(length(stock_names),
                                            missing_sample_events,
                                            simplify = FALSE)) %>%
    mutate(regName = rep(stock_names, each = nrow(missing_sample_events))) %>%
    select(id:season, regName, pres, catchReg)
  
  # combine and check for no zeros
  temp <- comp_trim %>%
    rbind(., infill_data)
  tab_in <- table(comp_trim$regName, comp_trim$month, comp_trim$catchReg)
  tab_out <- table(temp$regName, temp$month, temp$catchReg)
  
  # spread
  gsi_out <- temp %>% 
    mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
    pivot_wider(., names_from = regName, values_from = pres) %>%
    mutate_if(is.numeric, ~replace_na(., 0)) %>%
    droplevels() 
  
  ifelse(check_tables == FALSE, 
         return(gsi_out), 
         return(list(data = gsi_out, tables = list(original = tab_in, 
                                                   infilled = tab_out)))
         )
}


# Prep Data for TMB ------------------------------------------------------------
tmb_dat <- function(catch, comp_wide, fac_dat, fac_key, mod = "nb") {
  # fixed effects model matrices that either include or exclude effort
  sp_t_ch <- paste("~ catchReg + month")
  sp_t_form <- formula(sp_t_ch)
  sp_t_eff_ch <- paste("~ catchReg + month + eff_z + eff_z2")
  sp_t_eff_form <- formula(sp_t_eff_ch)
  
  fix_mm_gsi <- model.matrix(sp_t_form, data = comp_wide)
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

  # predictive model matrix for abundance based on fac key with the potential for
  # the addition of covariate effects
  mm_pred <- model.matrix(sp_t_eff_form, data = fac_key)

  if(!identical(colnames(mm_pred), colnames(fix_mm_c))) {
    warning("Mismatch between estimated and predicted abundance model matrices.")
  }
  
  # observed stock composition
  y_obs <- comp_wide %>% 
    select(-c(id:dummy_id)) %>% 
    as.matrix()
  head(y_obs)
  
  # other parameters
  n_groups <- ncol(y_obs)
  b1_n <- length(coef(m1))
  b2_n <- ncol(mm_pred) * (n_groups - 1) #each par est for each non-ref level
  # vectors of random effects
  fac1k <- fct_to_tmb_num(catch$year)
  fac2k <- fct_to_tmb_num(comp_wide$year)
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
 
  list(data = data, parameters = parameters)
}

# Fit Model --------------------------------------------------------------------

source(here::here("R", "functions", "clean_composition_dat.R"))
comp_wide <- clean_gsi(comp, month_range = month_range)

# make fixed effects factor key based on stock composition data
fac_dat <- comp_wide %>% 
  mutate(facs = as.factor(paste(as.character(catchReg), 
                                as.character(month), sep = "_")),
         facs = fct_reorder2(facs, 
                             desc(month_n),
                             desc(catchReg)),
         facs = fct_relevel(facs, "NWVI_10", after = 9),
         facs = fct_relevel(facs, "SWVI_10", after = Inf),
         facs_n = as.numeric(facs) - 1) %>% 
  select(catchReg, month, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)
fac_key_eff <- fac_key %>%
  mutate(eff_z = mean(catch$eff_z),
         eff_z2 = mean(catch$eff_z2))

source(here::here("R", "functions", "prep_tmb_dat.R"))
inputs <- tmb_dat(catch, comp_wide, fac_dat, fac_key_eff, mod = "nb")

## Make a function object
# compile(here::here("R", "southCoast", "tmb", "tweedie_multinomial_1re.cpp"))
# dyn.load(dynlib(here::here("R", "southCoast", "tmb", 
#                            "tweedie_multinomial_1re")))
# obj <- MakeADFun(data, parameters, random = c("z1_k", "z2_k"), 
#                  DLL = "tweedie_multinomial_1re")

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

saveRDS(ssdr, here::here("generatedData", "model_fits", "nbmult_ssdr_agg.RDS"))
# saveRDS(ssdr, here::here("generatedData", "model_fits", "twmult_ssdr_frB.RDS"))

# PREDICTIONS ------------------------------------------------------------------

#frB versions have a different subset of stocks and months at finer spatial scale
# ssdr <- readRDS(here::here("generatedData", "model_fits", "twmult_ssdr_frB.RDS"))
# ssdr <- readRDS(here::here("generatedData", "model_fits", "twmult_ssdr_agg.RDS"))

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
  # rename(stock = regName) %>% 
  filter(!is.na(stock))


# combined estimates of stock-specific CPUE
# note that effort is standardized differently between raw data and predictions,
# not apples to apples comparison
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

# estimates of aggregate CPUE (replace with simulated approach below)
# agg_abund <- data.frame(
#   facs_n = fac_key$facs_n,
#   raw_abund_est = log_pred[ , "Estimate"],
#   raw_abund_se = log_pred[ , "Std. Error"]) %>% 
#   mutate(
#     raw_mu = exp(raw_abund_est),
#     raw_abund_low = exp(raw_abund_est + (qnorm(0.025) * raw_abund_se)),
#     raw_abund_up = exp(raw_abund_est + (qnorm(0.975) * raw_abund_se))
#   ) %>% 
#   left_join(pred_ci, ., by = "facs_n")

# catch <- catch %>% 
#   group_by(month, catchReg) %>% 
#   mutate(month_eff = mean(boatDays)) %>% 
#   ungroup() %>% 
#   mutate(catch_z = catch / month_eff)
# 
# ggplot() +
#   geom_boxplot(data = catch %>% filter(!catch == 0), aes(x = month, y = cpue),  
#              alpha = 0.4) +
#   geom_pointrange(data = agg_abund, aes(x =  month, y = raw_mu,
#                                   ymin = raw_abund_low,
#                                   ymax = raw_abund_up), color= "red") +
#   facet_wrap(~ catchReg, nrow = 2, scales = "free_y") +
#   ggsidekick::theme_sleek()


## generate better comparison of abundance model predictions by using betas
## to calculate catch across range of effort values

#betas for abundance model
abund_b <- ssdr[rownames(ssdr) %in% "b1_j", ]
pred_catch <- fac_key %>% 
  mutate(
    #add model estimates
    int = abund_b[1],
    reg_b = case_when(
      catchReg == "NWVI" ~ 0,
      catchReg == "SWVI" ~ abund_b[2]
    ),
    month_b = rep(c(0, abund_b[3:11]), times = 2),
    eff_b = abund_b[12],
    eff2_b = abund_b[13]) %>% 
  split(., .$facs_n) %>% 
  map(., function(x) {
    mm <- catch %>% 
      filter(month == x$month) 
    expand_grid(x, z_eff = sample(mm$eff_z, size = 50, replace = T))
  }) %>% 
  bind_rows() %>% 
  #add estimates
  mutate(z_eff2 = z_eff^2,
         log_catch = int + reg_b + month_b + (eff_b * z_eff) + 
           (eff2_b * z_eff2),
         catch = exp(log_catch),
         dataset = "pred")

catch %>% 
  mutate(dataset = "obs") %>% 
  select(catch, dataset, catchReg, month) %>% 
  rbind(., 
        pred_catch %>% 
          select(catch, dataset, catchReg, month)
        ) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = month, y = catch, fill = dataset)) +
  facet_wrap(~ catchReg, nrow = 2, scales = "free_y") +
  ggsidekick::theme_sleek()
# looks good! some deviations because effort isn't stratified by region and 
# July data are wonky, but otherwise solid

## export plotting data for Rmd (version 1)
# list(catch = catch, pred_ci = pred_ci, raw_prop = raw_prop, raw_abund = raw_abund, 
#      raw_agg_abund = agg_abund) %>% 
#   saveRDS(., here::here("generatedData", "model_fits", "frB_plot_list.RDS"))

## export plotting data for Rmd (version 1)
list(catch = catch, pred_ci = pred_ci, raw_prop = raw_prop, raw_abund = raw_abund,
     pred_catch = pred_catch) %>%
  saveRDS(., here::here("generatedData", "model_fits", "reg3_plot_list.RDS"))


## estimates of effort effects on catch
n_betas <- length(coef(m1))
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

