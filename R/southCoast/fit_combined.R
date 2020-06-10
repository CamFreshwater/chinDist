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
    mutate(region_c = as.character(region),
           region = factor(region_c),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff)),
           eff_z2 = eff_z^2) %>% 
    arrange(region, month) 
}

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  rename(eff = boatDays, region = catchReg) %>% 
  clean_catch(.) %>% 
  # drop month where catch data are missing even though gsi samples available
  mutate(temp_strata = paste(month_n, region, sep = "_")) %>% 
  filter(!temp_strata == "7_SWVI") 

#recreational catch data
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "monthlyCatch_rec.RDS")) %>% 
  #group to get rid of adipose and released legal categories
  group_by(month, month_n, year, area, subarea, region, legal) %>% 
  summarize(catch = sum(mu_catch),
            eff = mean(mu_boat_trips),
            cpue = catch / eff) %>% 
  ungroup() %>% 
  mutate(temp_strata = paste(month_n, region, sep = "_"),
         region = abbreviate(region, minlength = 4))  %>% 
  clean_catch(.) %>%  
  #focus on legal fish 
  filter(legal == "legal")

tt <- rec_catch %>% 
  filter(legal == "legal") %>% 
  group_by(subarea, year, month) %>% 
  summarize(n_eff = length(unique(mu_boat_trips)))



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
  # mutate(temp_strata = paste(month, region, sep = "_")) %>% 
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
    mutate(region_c = as.character(region),
           region = as.factor(region))
  
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


## PREP MODEL INPUTS -----------------------------------------------------------

# function to infill genetics data and subset data based on months present
subset_comp <- function(x, threshold = 50) {
  # original structure only remove strata 
  # (SWITCH BACK ONCE PARS CAN BE FIXED IN TMB)
  # x %>% 
  #   mutate(temp_strata = paste(month_n, region, sep = "_")) %>% 
  #   group_by(temp_strata) %>% 
  #   mutate(nn = length(unique(id))) %>%
  #   filter(!nn < threshold) %>% 
  #   ungroup() %>% 
  #   droplevels() %>% 
  #   select(id, region, area, year, month, season, agg, agg_prob, pres, 
  #          temp_strata)
  
  # more conservative structure removes all months associated w/ missing strata
  empty_grid <- expand_grid(month = unique(x$month), 
                            region = unique(x$region)) %>% 
    mutate(id = "dum")
  dum <- x %>% 
    select(month, region, id) %>% 
    rbind(empty_grid) %>% 
    group_by(month, region) %>% 
    summarize(nn = length(unique(id))) %>%
    filter(nn < threshold) %>% 
    ungroup()
  
  x %>% 
    filter(!month %in% unique(dum$month)) %>%  
    droplevels() %>% 
    select(id, region, area, year, month, season, agg, agg_prob, pres, 
           temp_strata) 
}

# function to infill then widen composition data
infill_comp_dat <- function(x) {
  #make wide version of gsi data and infill 0 observations
  #has to occur for each region separately given differences in sample effort
  dum <- x %>% 
    group_by(region) %>% 
    nest() %>%
    transmute(dum_data = map(data, function(x) {
      expand.grid(month = unique(x$month),
                  agg = unique(x$agg),
                  pres = 1)
    })
    ) %>% 
    unnest(cols = c(dum_data))
  # add random years 
  rand_yrs <- sample(unique(x$year), size = nrow(dum), replace = TRUE)
  dum$year <- rand_yrs
  
  suppressWarnings(
    x %>%
      full_join(., dum, by = c("year", "month", "agg", "pres", "region")) %>%
      mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
      droplevels() 
  )
}

widen_comp_dat <- function(long_comp) {
  long_comp %>% 
    pivot_wider(., names_from = agg, values_from = pres) %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
}

# function to make predictive dataframes 
make_pred_dat <- function(dat) {
  dat %>%
    select(region, month) %>%
    distinct() %>%
    arrange(region, month)
}

# function to generate TMB data inputs from wide dataframes
gen_tmb_dat <- function(catch_dat, comp_dat, effort) {
  ## catch data
  yr_vec_catch <- as.numeric(as.factor(as.character(catch_dat$year))) - 1
  
  # model matrix for fixed effects
  fix_mm_catch <- model.matrix(catch_formula, catch_dat)
  
  # use average effort by month and area for predictions
  if (effort == "month") {
    mean_eff <- catch_dat %>% 
      group_by(region, month) %>%
      summarize(eff_z = mean(eff_z),
                eff_z2 = mean(eff_z2)) %>% 
      ungroup()
    pred_dat_catch <- make_pred_dat(catch_dat) %>% 
      left_join(., mean_eff, by = c("region", "month"))
    pred_mm_catch <- model.matrix(catch_formula, pred_dat_catch) 
  }  
  # set effort to zero for all predictions
  if (effort == "overall") {
    pred_dat_catch <- make_pred_dat(catch_dat)
    pred_mm_catch <- model.matrix(comp_formula, pred_dat_catch) %>%
      cbind(.,
            eff_z = rep(0, n = nrow(.)),
            eff_z2 = rep(0, n = nrow(.)))
  }
  
  
  ## composition data
  obs_comp <- comp_dat %>% 
    select(-c(id:dummy_id)) %>% 
    as.matrix()
  yr_vec_comp <- as.numeric(comp_dat$year) - 1
  fix_mm_comp <- model.matrix(comp_formula, comp_dat)
  
  pred_dat_comp <- make_pred_dat(comp_dat)
  pred_mm_comp <- model.matrix(~ region + month, pred_dat_comp)
  
  list(
    #abundance input data
    y1_i = catch_dat$catch,
    X1_ij = fix_mm_catch,
    factor1k_i = yr_vec_catch,
    nk1 = length(unique(yr_vec_catch)),
    X1_pred_ij = pred_mm_catch,
    #composition input data
    y2_ig = obs_comp,
    X2_ij = fix_mm_comp,
    factor2k_i = yr_vec_comp,
    nk2 = length(unique(yr_vec_comp)),
    X2_pred_ij = pred_mm_comp
  )
}

# function to generate TMB parameters
gen_tmb_par <- function(catch_dat, tmb_data) {
  #fit simple model to get initial values for abundance pars
  full_catch_formula <- as.formula(paste("log(catch + 0.0001)", "~", 
                                         as.character(catch_formula)[2],
                   sep = " "))
  m1 <- lm(full_catch_formula, data = catch_dat)
  
  list(
    #abundance parameters
    b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
    log_phi = log(1.5),
    z1_k = rep(0, tmb_data$nk1),
    log_sigma_zk1 = log(0.25),
    #composition parameters
    b2_jg = matrix(0, 
                   nrow = ncol(tmb_data$X2_pred_ij), 
                   ncol = ncol(tmb_data$y2_ig) - 1),
    z2_k = rep(0, times = tmb_data$nk2),
    log_sigma_zk2 = log(0.25)
  )
}

# relevant formulas
comp_formula <- as.formula("~ region + month")
catch_formula <- as.formula("~ region + month + eff_z + eff_z2")

# join all data into a tibble then generate model inputs
dat <- tibble(
  dataset = c("pst_rec", "pst_comm", "can_rec", "can_comm"),
  catch = list(rec_catch, comm_catch, rec_catch, comm_catch),
  comp = list(rec_pst_comp, comm_pst_comp, rec_can_comp, comm_can_comp)
  ) %>% 
  transmute(
    dataset,
    #subset composition data based on number of samples collected
    long_comp = map(comp, subset_comp, threshold = 50),
    # subset catch based on composition data
    catch = map2(catch, long_comp, function (x, y) {
      x %>% 
        filter(temp_strata %in% y$temp_strata) %>% 
        droplevels()
    }),
    # infill then widen composition data
    infill_comp = map(long_comp, infill_comp_dat),
    comp = map(infill_comp, widen_comp_dat),
    # predictive dataframes for catch
    pred_catch_dat = map(catch, make_pred_dat),
    # predictive dataframes for comp
    pred_comp_dat = map(comp, make_pred_dat),
    #add tmb data assuming month-area specific effort predictions
    # tmb_data = map2(catch, comp, gen_tmb_dat, effort = "month"),
    #add tmb data assuming average effort predictions
    tmb_data = map2(catch, comp, gen_tmb_dat, effort = "overall"),
    #add tmb parameters
    tmb_pars = map2(catch, tmb_data, gen_tmb_par)
  ) 

# check frequency tables of composition data
map(dat$infill_comp, function (x) {
  table(x$agg, x$month, x$region)
})


## EXAMINE RAW DATA ------------------------------------------------------------

catch_list <- list(dat$catch[[1]], dat$catch[[2]])
map(catch_list, function(x) {
    cpue <- ggplot(x) +
      geom_boxplot(aes(x = month, y = cpue)) +
      facet_wrap(~region) +
      ggsidekick::theme_sleek()
    effort <- x %>% 
      select(region, month, eff) %>% 
      distinct() %>% 
      ggplot(.) +
      geom_boxplot(aes(x = month, y = eff)) +
      facet_wrap(~region) +
      ggsidekick::theme_sleek()
    cowplot::plot_grid(cpue, effort, nrow = 2)
  })


## FIT MODELS ------------------------------------------------------------------

compile(here::here("src", "nb_multinomial_1re_v2.cpp"))
dyn.load(dynlib(here::here("src", "nb_multinomial_1re_v2")))

fit_mod <- function(tmb_data, tmb_pars) {
  # tmb_map <- list(b2_jg = factor(tmb_pars$b2_jg))
  obj <- MakeADFun(tmb_data, tmb_pars, random = c("z1_k", "z2_k"), 
                   DLL = "nb_multinomial_1re_v2"
                   #, map = tmb_map
                   )
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  sdreport(obj)
  # ssdr <- summary(sdr)
}

dat2 <- dat %>% 
  mutate(sdr = map2(tmb_data, tmb_pars, fit_mod),
         ssdr = map(sdr, summary))
saveRDS(dat2, here::here("generated_data", "model_fits", "combined_model.RDS"))


## GENERATE PREDICTIONS --------------------------------------------------------

## Abundance predictions
abund_pred <- ssdr[rownames(ssdr) %in% "pred_abund", ] #log pred of abundance
abund_b <- ssdr[rownames(ssdr) %in% "b1_j", ]
catch <- dat$catch[[4]]
pred_catch <- dat$pred_catch_dat[[4]]

# plot predicted abundance with average effort (i.e. set to 0)
pred_ci <- data.frame(pred_est = abund_pred[ , "Estimate"],
                      pred_se =  abund_pred[ , "Std. Error"]) %>%
  cbind(pred_catch, .) %>% 
  mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
         pred_up = pred_est + (qnorm(0.975) * pred_se))
  
real_preds <- ggplot() +
  geom_pointrange(data = pred_ci, aes(x = as.factor(month), y = pred_est,
                                      ymin = pred_low, ymax = pred_up)) +
  labs(x = "month", y = "predicted real catch (mean effort)") +
  ggsidekick::theme_sleek() +
  facet_wrap(~region)

# plot predicted vs observed abundance accounting for month-area specific effort
pred_eff <- catch %>% 
  select(region, month, eff_z, eff_z2) %>% 
  group_by(region, month) %>% 
  sample_n(., size = 50, replace = TRUE) %>% 
  ungroup()
pred_mm2 <- model.matrix(catch_formula, pred_eff)
pred_catch2 <- pred_mm2 %*% abund_b 
pred_plot <- pred_eff %>% 
  mutate(log_catch = pred_catch2[ , 'Estimate'],
         catch = exp(log_catch),
         dataset = "pred")

var_effort_preds <- catch %>% 
  mutate(dataset = "obs") %>% 
  select(catch, dataset, region, month) %>% 
  rbind(., 
        pred_plot %>% 
          select(catch, dataset, region, month)
  ) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = month, y = catch, fill = dataset)) +
  facet_wrap(~ region, nrow = 2, scales = "free_y") +
  ggsidekick::theme_sleek() +
  labs(fill = "Data")


## Composition predictions
comp_pred <- ssdr[rownames(ssdr) %in% "pred_probs", ]
comp_abund_pred <- ssdr[rownames(ssdr) %in% "pred_abund_mg", ]
comp <- dat$comp[[4]]
comp_long <- dat$long_comp[[4]]
y_obs <- dat$tmb_data[[4]]$y2_ig
stk_names <- colnames(y_obs)
N <- nrow(y_obs)
k <- ncol(y_obs)

pred_comp_dat <- purrr::map_dfr(seq_len(k), ~ dat$pred_comp_dat[[4]])
pred_comp_ci <- data.frame(stock = as.character(rep(stk_names, 
                                               each = nrow(dat$pred_comp_dat[[4]]))),
                      pred_prob_est = comp_pred[ , "Estimate"],
                      pred_prob_se =  comp_pred[ , "Std. Error"]) %>% 
  cbind(pred_comp_dat, .) %>%
  mutate(pred_prob_low = pred_prob_est + (qnorm(0.025) * pred_prob_se),
         pred_prob_up = pred_prob_est + (qnorm(0.975) * pred_prob_se),
         comp_abund_est = comp_abund_pred[ , "Estimate"],
         comp_abund_se =  comp_abund_pred[ , "Std. Error"],
         comp_abund_low = comp_abund_est + (qnorm(0.025) * comp_abund_se),
         comp_abund_up = comp_abund_est + (qnorm(0.975) * comp_abund_se)) 

# calculate raw proportion data for comparison
raw_prop <- comp_long %>% 
  group_by(region, month, year, agg) %>%
  summarize(samp_g = length(unique(id))) %>% 
  group_by(region, month, year) %>%
  mutate(samp_total = sum(samp_g)) %>% 
  ungroup() %>% 
  mutate(samp_g_ppn = samp_g / samp_total,
         stock = fct_reorder(agg, desc(samp_g_ppn))) 

comp_preds <- ggplot() +
  geom_point(data = raw_prop,
             aes(x = month, y = samp_g_ppn, fill = region),
             shape = 21, alpha = 0.4, position = position_dodge(0.6)) +
  geom_pointrange(data = pred_comp_ci,
                  aes(x = month, y = pred_prob_est,
                      ymin = pred_prob_low, ymax = pred_prob_up,
                      fill = region),
                  shape = 21, size = 0.4, position = position_dodge(0.6)) +
  labs(y = "Probability", x = "Month") +
  facet_wrap(~ stock) +
  ggsidekick::theme_sleek()


## Stock-specific abundance predictions
raw_abund <- catch %>% 
  group_by(region, month, year) %>%
  summarize(sum_catch = sum(catch),
            sum_effort = sum(eff),
            agg_cpue = sum_catch / sum_effort) %>% 
  ungroup() %>% 
  mutate(month = as.character(month)) %>% 
  left_join(., raw_prop, by = c("region", "month", "year")) %>% 
  mutate(catch_g = samp_g_ppn * sum_catch,
         cpue_g = samp_g_ppn * agg_cpue,
         region = as.factor(region),
         month = fct_relevel(as.factor(month), "10", after = Inf)
         ) %>% 
  filter(!is.na(stock))


# combined estimates of stock-specific CPUE
# note that effort is standardized differently between raw data and predictions,
# not apples to apples comparison
ggplot() +
  geom_point(data = raw_abund, aes(x = month, y = cpue_g, fill = region),
             shape = 21, alpha = 0.3, position = position_dodge(0.6)) +
  geom_pointrange(data = pred_comp_ci,
                  aes(x = month, y = comp_abund_est, ymin = comp_abund_low,
                      ymax = comp_abund_up, fill = region),
                  shape = 21, position = position_dodge(0.6)) +
  facet_wrap(~stock, ncol = 2, scales = "free_y") +
  scale_fill_viridis_d(option = "C", name = "Catch Region") +
  labs(x = "Month", y = "Predicted Catch") +
  ggsidekick::theme_sleek() +
  theme(legend.position="top")

ggplot() +
  geom_point(data = raw_abund, aes(x = month, y = cpue_g, fill = region),
             shape = 21, alpha = 0.7, position = position_dodge(0.6)) +
  facet_wrap(~stock, ncol = 2, scales = "free_y") +
  scale_fill_viridis_d(option = "C", name = "Catch Region") +
  labs(x = "Month", y = "Predicted Catch") +
  ggsidekick::theme_sleek() +
  theme(legend.position="top")




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
