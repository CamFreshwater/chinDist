## Combined model fits
# April 3, 2020
# Fit combined multionomial/tweedie or multinomial/nb model to 
# stock composition and abundance data
# aggregate probability reflects summed probabilities of a given region of 
# origin for a given individual

library(tidyverse)
library(TMB)
library(ggplot2)


# CLEAN CATCH -----------------------------------------------------------------

# clean function 
clean_catch <- function(dat) {
  dat %>% 
    filter(!is.na(eff),
           !eff == "0") %>% 
    mutate(reg_f = factor(region),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff)),
           eff_z2 = eff_z^2,
           eff_z3 = eff_z^3) %>% 
    arrange(reg_f, month) 
}

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  rename(eff = boatDays, region = catchReg) %>% 
  clean_catch(.) %>% 
  # drop month where catch data are missing even though gsi samples available
  mutate(temp_strata = paste(month, region, sep = "_")) %>% 
  filter(!temp_strata == "7_SWVI") 

#recreational catch data
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "monthlyCatch_rec.RDS")) %>% 
  # drop months that lack data from all regions
  rename(catch = mu_catch, eff = mu_boat_trips) %>% 
  mutate(cpue = catch / eff,
         region = abbreviate(region, minlength = 4)) %>% 
  clean_catch(.) %>%  
  #focus on legal fish 
  filter(legal == "legal")


# CLEAN GENETICS  --------------------------------------------------------------

# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  #focus only on legal sized fish
  filter(legal == "legal") %>% 
  mutate(region = abbreviate(region, minlength = 4))
# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds")) %>% 
  # drop month where catch data are missing
  mutate(temp_strata = paste(month, region, sep = "_")) %>% 
  filter(!temp_strata == "7_SWVI")


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
    select(-max_assignment) %>% 
    mutate(reg_f = as.factor(region))
  
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
      ),
      reg1 = case_when(
        Region1Name %in% c("Fraser_Fall", "ECVI", "Fraser_Summer_4.1", 
                           "Fraser_Spring_4.2", "Fraser_Spring_5.2", 
                           "Fraser_Summer_5.2", "WCVI", "SOMN") ~ Region1Name,
        TRUE ~ "Other"
      )
    )
}

## combine genetics and catch data
#large scale PST groups 
rec_pst_comp <- rec %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., rec, thresh = 0.75)
comm_pst_comp <- comm %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., comm, thresh = 0.75) 
#aggregate with Canadian CU focus
rec_can_comp <- rec %>% 
  pool_aggs() %>% 
  group_by(id, reg1) %>%
  calc_max_prob(., rec, thresh = 0.75) 
comm_can_comp <- comm %>% 
  pool_aggs() %>% 
  group_by(id, reg1) %>%
  calc_max_prob(., comm, thresh = 0.75) 

dat <- tibble(
  dataset = c("pst_rec", "pst_comm", "can_rec", "can_comm"),
  catch = list(rec_catch, comm_catch, rec_catch, comm_catch),
  comp = list(rec_pst_comp, comm_pst_comp, rec_can_comp, comm_can_comp)
)

dat_list <- list("pst_rec" = rec_pst, "pst_comm" = comm_pst, 
                 "can_rec" = rec_can, "can_comm" = comm_can) %>% 
  # remove months where too few GSI samples collected
  map(., function (x, threshold = 50) {
    comp_out <- x$comp %>% 
      mutate(temp_strata = paste(month_n, region, sep = "_")) %>% 
      group_by(temp_strata) %>% 
      mutate(nn = length(unique(id))) %>%
      filter(!nn < threshold) %>% 
      ungroup() %>% 
      droplevels() %>% 
      select(id, reg_f, area, year, month, season, agg, agg_prob, pres, 
             temp_strata)
    catch_out <- x$catch %>% 
      mutate(temp_strata = paste(month, region, sep = "_")) %>% 
      filter(temp_strata %in% comp_out$temp_strata) %>% 
      droplevels()
    
    #make wide version of gsi data and infill 0 observations
    #has to occur for each region separately given differences in sample effort
    comp_out %>% 
      split(., .$reg_f) %>% 
      map(., function(x) {
        expand.grid(month = unique(x$month),
                    agg = unique(x$agg),
                    pres = 1)
      }) %>% 
    
    dum <- expand.grid(
      month = unique(comp_out$month),
      reg_f = unique(comp_out$reg_f),
      agg = unique(comp_out$agg),
      pres = 1)
    rand_yrs <- sample(unique(comp_out$year), size = nrow(dum), replace = TRUE)
    dum$year <- rand_yrs
    
    gsi_trim_no0 <- comp_out %>%
      full_join(., dum, by = c("year", "month", "agg", "pres", "reg_f")) %>%
      mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
      droplevels() 
    freq_table <- table(gsi_trim_no0$agg, gsi_trim_no0$month, gsi_trim_no0$reg_f)
    
    comp_wide <- gsi_trim_no0 %>% 
      pivot_wider(., names_from = agg, values_from = pres) %>% 
      mutate_if(is.numeric, ~replace_na(., 0))
    
    list("catch" = catch_out, "comp" = comp_wide "comp_long" =  comp_out)
  })


## GENERATE INPUTS -------------------------------------------------------------

list_in$comp %>% 
  group_by(reg_f) %>% 
  tally()
list_in$catch %>% 
  group_by(reg_f) %>% 
  tally()

comp_formula <- paste("~ catchReg + month")
catch_formula <- paste("~ catchReg + month + eff_z + eff_z2")

prep_data <- function(catch_formula, comp_formula, list_in) {
  
}

prep_data <- function(catch, data_type = NULL) {
  
  yr_vec <- as.numeric(as.factor(as.character(catch$year))) - 1
  
  # model matrix for fixed effects
  fix_mm <- model.matrix(~ reg_f + month + eff_z + eff_z2, catch)
  
  # Average effort for predictions
  # pred_eff <- catch %>% 
  #   mutate(facs = as.character(paste(reg_f, month_n, sep = "_"))) %>% 
  #   group_by(facs) %>% 
  #   summarize(eff_z = mean(eff_z),
  #             eff_z2 = mean(eff_z2))
  
  # Factor key of unique combinations to generate predictions
  fac_key <- catch %>%
    select(reg_f, month_n) %>%
    distinct() %>%
    mutate(month = as.factor(month_n),
           facs = paste(reg_f, month_n, sep = "_"),
           facs = fct_reorder2(facs, reg_f, 
                               desc(month_n)),
           facs_n = as.numeric(as.factor(facs)) - 1
    ) %>%
    arrange(facs_n) %>% 
    left_join(., pred_eff, by = "facs") %>% 
    mutate(facs = as.factor(facs))
  
  #mm_pred <- model.matrix(~ reg_f + month + eff_z + eff_z2, fac_key)
  mm_pred <- model.matrix(~ reg_f + month + eff_z + eff_z2, fac_key)%>% 
    cbind(.,
          eff_z = rep(0, n = nrow(.)),
          eff_z2 = rep(0, n = nrow(.)))
  
  data <- list(y1_i = catch$catch,
               X1_ij = fix_mm,
               factor1k_i = yr_vec,
               nk1 = length(unique(yr_vec)),
               X1_pred_ij = mm_pred
  )
  
  # Fit simple model to initialize tmb 
  m1 <- lm(log(catch + 0.0001) ~ reg_f + month + eff_z + eff_z2, data = catch)
  
  parameters = list(
    b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
    log_phi = log(1.5),
    z1_k = rep(0, length(unique(yr_vec))),
    log_sigma_zk1 = log(0.25)
  )
  
  if (is.null(data_type)) {
    data_type <- unique(catch$legal)
  }
  
  list("fix_mm" = fix_mm, "fac_key" = fac_key, "mm_pred" = mm_pred, 
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "input_data" = catch)
}

prep_gsi <- function(gsi_in, month_range = c(1, 12), data_type) {
  gsi_trim <- gsi_in %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(id, region, area, year, month, season, agg, agg_prob, pres)
  
  # dummy dataset to replace missing values 
  dum <- expand.grid(
    month = unique(gsi_trim$month),
    region = unique(gsi_trim$region),
    agg = unique(gsi_trim$agg),
    pres = 1)
  rand_yrs <- sample(unique(gsi_trim$year), size = nrow(dum), replace = TRUE)
  dum$year <- rand_yrs
  
  gsi_trim_no0 <- gsi_trim %>%
    full_join(., dum, by = c("year", "month", "agg", "pres", "region")) %>%
    mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
    droplevels() 
  
  gsi_wide <- gsi_trim_no0 %>% 
    pivot_wider(., names_from = agg, values_from = pres) %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
  raw_freq_table <- table(gsi_trim$agg, gsi_trim$month, gsi_trim$region)
  freq_table <- table(gsi_trim_no0$agg, gsi_trim_no0$month, gsi_trim_no0$region)
  
  y_obs <- gsi_wide %>% 
    select(-c(id:dummy_id)) %>% 
    as.matrix()
  head(y_obs)
  
  yr_vec <- as.numeric(gsi_wide$year) - 1
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
  
  list("fix_mm" = fix_mm, #"pred_dat" = pred_dat, 
       "fac_key" = fac_key,
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "long_data" = gsi_trim, "long_data_no0" = gsi_trim_no0, 
       "wide_data" = gsi_wide, "freq_table" = freq_table, 
       "raw_freq_table" = raw_freq_table)
}






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
inputs <- tmb_dat(catch, comp_wide, fac_dat, fac_key_eff, mod = "nb")
#check predictive matrix
head(inputs$data$X1_ij)
head(inputs$data$X1_pred_ij)


# Fit Model --------------------------------------------------------------------

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

# saveRDS(ssdr, here::here("generatedData", "model_fits", "nbmult_ssdr_pst.RDS"))
# saveRDS(ssdr, here::here("generatedData", "model_fits", "nbmult_ssdr_frB.RDS"))


# PREDICTIONS ------------------------------------------------------------------

#frB versions have a different subset of stocks and months at finer spatial scale
ssdr <- readRDS(here::here("generatedData", "model_fits", "nbmult_ssdr_pst.RDS"))
# ssdr <- readRDS(here::here("generatedData", "model_fits", "twmult_ssdr_agg.RDS"))

log_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] #log pred of abundance
logit_probs <- ssdr[rownames(ssdr) %in% "logit_pred_prob", ] #logit probs of each category
pred_abund <- ssdr[rownames(ssdr) %in% "pred_abund_mg", ] #pred abundance of each category

comp_trim <- comp_wide_l$long_data
n_groups <- length(unique(comp_trim$regName))

pred_ci <- data.frame(stock = as.character(rep(unique(comp_trim$regName),
                                               each = 
                                                 length(unique(fac_key_eff$facs_n)))), 
                      logit_prob_est = logit_probs[ , "Estimate"],
                      logit_prob_se =  logit_probs[ , "Std. Error"]) %>%
  mutate(facs_n = rep(fac_key_eff$facs_n, times = n_groups)) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_low = plogis(logit_prob_est +
                                  (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est +
                                 (qnorm(0.975) * logit_prob_se)),
         abund_est = pred_abund[ , "Estimate"],
         abund_se =  pred_abund[ , "Std. Error"],
         abund_low = abund_est + (qnorm(0.025) * abund_se),
         abund_up = abund_est + (qnorm(0.975) * abund_se)) %>%
  left_join(., fac_key_eff, by = c("facs_n")) %>% 
  mutate(stock = fct_reorder(stock, desc(pred_prob)))

# calculate raw summary data for comparison
raw_prop <- comp %>% 
  filter(month_n %in% month_range) %>% 
  group_by(catchReg, month, year, regName) %>%
  summarize(samp_g = length(unique(id))) %>% 
  group_by(catchReg, month, year) %>%
  mutate(samp_total = sum(samp_g)) %>% 
  ungroup() %>% 
  mutate(samp_g_ppn = samp_g / samp_total,
         stock = fct_reorder(regName, desc(samp_g_ppn))) 

raw_abund <- catch %>% 
  group_by(catchReg, month, year) %>%
  summarize(sum_catch = sum(catch),
            sum_effort = sum(boatDays),
            agg_cpue = sum_catch / sum_effort) %>% 
  ungroup() %>% 
  mutate(month = as.character(month)) %>% 
  left_join(., raw_prop, by = c("catchReg", "month", "year")) %>% 
  mutate(catch_g = samp_g_ppn * sum_catch,
         cpue_g = samp_g_ppn * agg_cpue,
         catchReg = as.factor(catchReg), 
         month = fct_relevel(as.factor(month), "10", after = Inf)) %>% 
  filter(!is.na(stock))


# combined estimates of stock-specific CPUE
# note that effort is standardized differently between raw data and predictions,
# not apples to apples comparison
ggplot() +
  geom_point(data = raw_abund, aes(x = month, y = cpue_g, fill = catchReg),
             shape = 21, alpha = 0.3, position = position_dodge(0.6)) +
  geom_pointrange(data = pred_ci, 
                  aes(x = month, y = abund_est, ymin = abund_low, 
                      ymax = abund_up, fill = catchReg), 
                  shape = 21, position = position_dodge(0.6)) +
  facet_wrap(~stock, ncol = 2, scales = "free_y") +
  scale_fill_viridis_d(option = "C", name = "Catch Region") +
  labs(x = "Month", y = "Predicted Catch") +
  ggsidekick::theme_sleek() +
  theme(legend.position="top")

# estimates of stock compostion
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

# estimates of aggregate CPUE (replace with simulated approach below)
# agg_abund <- data.frame(
#   facs_n = fac_key_eff$facs_n,
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
pred_catch <- fac_key_eff %>% 
  mutate(
    #add model estimates
    int = abund_b[1],
    reg_b = case_when(
      catchReg == "NWVI" ~ 0,
      catchReg == "SWVI" ~ abund_b[2]
    ),
    month_b = rep(c(0, abund_b[3:(nrow(abund_b) - 2)]), times = 2),
    eff_b = abund_b[nrow(abund_b) - 1],
    eff2_b = abund_b[nrow(abund_b)]) %>% 
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
  ggsidekick::theme_sleek() +
  labs(fill = "Data")
# looks good! some deviations because effort isn't stratified by region and 
# July data are wonky, but otherwise solid

## export plotting data for Rmd 
list(catch = catch, pred_ci = pred_ci, raw_prop = raw_prop, raw_abund = raw_abund,
     pred_catch = pred_catch) %>%
  saveRDS(., here::here("generatedData", "model_fits", "pst_plot_list.RDS"))


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

# temp subset of stocks for comparison with cwt

stk_subset <- c("FR-early", "FR-late", "PSD", "CR-tule", "CR-bright", "WCVI")
plot_list <- map(list(pred_ci, raw_prop, raw_abund), function(x) 
  x %>% filter(stock %in% stk_subset)
)

prop_plot <- ggplot() +
  geom_pointrange(data = plot_list[[1]], aes(x = month, y = pred_prob, 
                                      ymin = pred_prob_low,
                                      ymax = pred_prob_up),
                  col = "red") +
  geom_point(data = plot_list[[2]],
             aes(x = month, y = samp_g_ppn),
             alpha = 0.4) +
  facet_wrap(stock ~ catchReg, nrow = n_groups, scales = "free_y") +
  ggsidekick::theme_sleek()

pdf(here::here("figs", "model_pred", "gsi_nb_prop_pred.pdf"))
prop_plot
dev.off()
