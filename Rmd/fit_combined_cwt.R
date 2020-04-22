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
  ) %>% 
  droplevels() %>% 
  mutate(eff_z = as.numeric(scale(boatDays)),
         eff_z2 = eff_z^2
  ) %>% 
  arrange(catchReg, area, month)


# Import Genetics --------------------------------------------------------------
# Function to trim and infill comp dataset
clean_comp <- function(comp, month_range, check_tables = FALSE) {
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
  comp_out <- temp %>% 
    mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
    pivot_wider(., names_from = regName, values_from = pres) %>%
    mutate_if(is.numeric, ~replace_na(., 0)) %>%
    droplevels() 
  
  ifelse(check_tables == FALSE, 
         return(comp_out), 
         return(list(data = comp_out, tables = list(original = tab_in, 
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

comp_wide <- clean_comp(comp, month_range = month_range)

# make fixed effects factor key based on stock composition data
fac_dat <- comp_wide %>% 
  mutate(facs = as.factor(paste(as.character(catchReg), 
                                as.character(month), sep = "_")),
         facs = fct_reorder2(facs, 
                             desc(month_n),
                             desc(catchReg)),
         facs = fct_relevel(facs, "NWVI_1", "NWVI_2", 
                            "NWVI_3", "NWVI_4", "NWVI_5", "NWVI_6",  "NWVI_7",
                            "NWVI_8", "NWVI_9", "NWVI_10", "NWVI_11", "NWVI_12", 
                            "SWVI_1",  "SWVI_2", "SWVI_3",  "SWVI_4", "SWVI_5", 
                            "SWVI_6",  "SWVI_7", 
                            "SWVI_8",  "SWVI_9",  "SWVI_10", "SWVI_11", "SWVI_12"),
         facs_n = as.numeric(facs) - 1) %>% 
  select(catchReg, month, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)
fac_key_eff <- fac_key %>%
  mutate(eff_z = mean(catch$eff_z),
         eff_z2 = mean(catch$eff_z2))

inputs <- tmb_dat(catch, comp_wide, fac_dat, fac_key_eff, mod = "nb")

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