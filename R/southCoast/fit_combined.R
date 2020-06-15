## Combined model fits
# April 3, 2020
# Fit combined multionomial/tweedie or multinomial/nb model to 
# stock composition and abundance data
# aggregate probability reflects summed probabilities of a given region of 
# origin for a given individual

library(tidyverse)
library(TMB)
library(ggplot2)
library(grid)
library(gridExtra)

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
  #group by subarea to get rid of adipose and released legal categories
  group_by(month, month_n, year, area, subarea, region, legal) %>% 
  mutate(subarea_catch = sum(mu_catch),
         subarea_eff = mean(mu_boat_trips)) %>% 
  #group by area to be consistent with commercial data
  group_by(month, month_n, year, area, region, legal) %>% 
  summarize(catch = sum(subarea_catch),
            eff = sum(subarea_eff),
            cpue = catch / eff) %>% 
  ungroup() %>% 
  mutate(temp_strata = paste(month_n, region, sep = "_"),
         region = abbreviate(region, minlength = 4))  %>% 
  clean_catch(.) %>%  
  #focus on legal fish 
  filter(legal == "legal")

# confirm one estimate of effort per month-subarea
# tt <- rec_catch %>%
#   filter(legal == "legal") %>%
#   group_by(subarea, year, month) %>%
#   summarize(n_eff = length(unique(mu_boat_trips)))
# range(tt$n_eff)


# CLEAN GENETICS  --------------------------------------------------------------

# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  #focus only on legal sized fish
  filter(legal == "legal") %>% 
  mutate(region_c = as.character(region),
         region = abbreviate(region, minlength = 4)) %>% 
  select(id:region, region_c, area, area_n, subarea:month_n, adj_prob:pst_agg)
# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds")) %>% 
  # drop month where catch data are missing
  filter(!temp_strata == "7_SWVI") %>% 
  mutate(region_c = as.character(region)) %>% 
  select(id:region, region_c, area, area_n, year:month_n, adj_prob:pst_agg)


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
                select(id:month_n) %>% 
                distinct(),
              by = "id") %>% 
    select(-max_assignment) %>% 
    mutate(region = as.factor(region))
  
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


## EXAMINE RAW DATA ------------------------------------------------------------

#catch and effort
catch_list <- list(comm_catch, rec_catch)
purrr::map(catch_list, function(x) {
  cpue <- ggplot() +
    geom_boxplot(data = x, aes(x = month, y = cpue)) +
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

#genetics samples collected
comp_list <- list(rec_pst_comp, comm_pst_comp)
purrr::map(comp_list, function(x) {
  x %>% 
    group_by(region, month, year) %>% 
    tally() %>% 
    ggplot(.) +
    geom_tile(aes(x = month, y = region, fill = n)) +
    facet_wrap(~year) +
    ggsidekick::theme_sleek()
})

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
    select(id, region, region_c, area, year, month, season, agg, agg_prob, pres, 
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

# function to calculate raw proportions based on observed composition data
make_raw_prop_dat <- function(comp_long) {
  comp_long %>% 
    group_by(region, month, year, agg) %>%
    summarize(samp_g = length(unique(id))) %>% 
    group_by(region, month, year) %>%
    mutate(samp_total = sum(samp_g)) %>% 
    ungroup() %>% 
    mutate(samp_g_ppn = samp_g / samp_total,
           stock = fct_reorder(agg, desc(samp_g_ppn))) 
}

# function to calculate raw stock-specific abundance observations
make_raw_abund_dat <- function(catch, raw_prop) {
  catch %>% 
    group_by(region, month, year) %>%
    summarize(sum_catch = sum(catch),
              sum_effort = sum(eff),
              agg_cpue = sum_catch / sum_effort) %>% 
    ungroup() %>% 
    mutate(month = as.character(month)) %>% 
    left_join(., raw_prop, by = c("region", "month", "year")) %>% 
    mutate(catch_g = samp_g_ppn * sum_catch,
           cpue_g = samp_g_ppn * agg_cpue,
           region = as.factor(region)) %>% 
    filter(!is.na(stock))
}

# function to make predictive dataframes 
make_pred_dat <- function(dat, catch = TRUE) {
  if (catch == TRUE) {
    dat %>%
      select(area, month, region) %>%
      distinct() %>%
      arrange(area, month, region)
  } else {
    dat %>%
      select(region, month) %>%
      distinct() %>%
      arrange(region, month)
  }
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
      group_by(area, month) %>%
      summarize(eff_z = mean(eff_z),
                eff_z2 = mean(eff_z2)) %>% 
      ungroup()
    pred_dat_catch <- make_pred_dat(catch_dat) %>% 
      left_join(., mean_eff, by = c("area", "month")) %>% 
      mutate(reg_month = paste(region, month, sep = "_"))
    pred_mm_catch <- model.matrix(catch_formula, pred_dat_catch) 
  }  
  # set effort to zero for all predictions
  if (effort == "overall") {
    pred_dat_catch <- make_pred_dat(catch_dat) %>% 
      mutate(reg_month = paste(region, month, sep = "_"))
    pred_mm_catch <- model.matrix(~ area + month, pred_dat_catch) %>%
      cbind(.,
            eff_z = rep(0, n = nrow(.)),
            eff_z2 = rep(0, n = nrow(.)))
  }
  
  # construct factor key for regional aggregates associated with areas
  grouping_vec <- as.numeric(as.factor(as.character(pred_dat_catch$reg_month))) - 1
  grouping_key <- unique(grouping_vec)
  
  
  ## composition data
  obs_comp <- comp_dat %>% 
    select(-c(id:dummy_id)) %>% 
    as.matrix()
  yr_vec_comp <- as.numeric(comp_dat$year) - 1
  fix_mm_comp <- model.matrix(comp_formula, comp_dat)
  
  pred_dat_comp <- make_pred_dat(comp_dat, catch = FALSE)
  pred_mm_comp <- model.matrix(comp_formula, pred_dat_comp)
  
  list(
    #abundance input data
    y1_i = catch_dat$catch,
    X1_ij = fix_mm_catch,
    factor1k_i = yr_vec_catch,
    nk1 = length(unique(yr_vec_catch)),
    X1_pred_ij = pred_mm_catch,
    pred_factor2k_h = grouping_vec,
    pred_factor2k_levels = grouping_key,
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
catch_formula <- as.formula("~ area + month + eff_z + eff_z2")

# join all data into a tibble then generate model inputs
dat <- tibble(
  dataset = c("pst_rec", "pst_comm", "can_rec", "can_comm"),
  catch = list(rec_catch, comm_catch, rec_catch, comm_catch),
  comp = list(rec_pst_comp, comm_pst_comp, rec_can_comp, comm_can_comp)
  ) %>% 
  transmute(
    dataset,
    #subset composition data based on number of samples collected
    long_comp = map(comp, subset_comp, threshold = 75),
    raw_prop = map(long_comp, gen_raw_prop),
    # subset catch based on composition data
    catch = map2(catch, long_comp, function (x, y) {
      x %>% 
        filter(temp_strata %in% y$temp_strata) %>% 
        droplevels()
    }),
    # raw estimates of composition and stock-specific abundance for plotting
    raw_prop = map(long_comp, make_raw_prop_dat),
    raw_abund = map2(catch, raw_prop, make_raw_abund_dat),
    # infill then widen composition data
    infill_comp = map(long_comp, infill_comp_dat),
    comp = map(infill_comp, widen_comp_dat),
    #add tmb data assuming average effort predictions
    tmb_data = map2(catch, comp, gen_tmb_dat, effort = "overall"),
    #add tmb parameters
    tmb_pars = map2(catch, tmb_data, gen_tmb_par),
    #add tmb data assuming month-area specific effort predictions
    #pars different due to different effort inputs
    tmb_data2 = map2(catch, comp, gen_tmb_dat, effort = "month"),
    tmb_pars2 = map2(catch, tmb_data2, gen_tmb_par)
    )


# check frequency tables of composition data
map(dat$infill_comp, function (x) {
  table(x$agg, x$month, x$region)
})


## FIT MODELS ------------------------------------------------------------------

compile(here::here("src", "nb_multinomial_1re_v3.cpp"))
dyn.load(dynlib(here::here("src", "nb_multinomial_1re_v3")))

fit_mod <- function(tmb_data, tmb_pars) {
  # tmb_map <- list(b2_jg = factor(tmb_pars$b2_jg))
  obj <- MakeADFun(tmb_data, tmb_pars, random = c("z1_k", "z2_k"), 
                   DLL = "nb_multinomial_1re_v3"
                   #, map = tmb_map
                   )
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdreport(obj)
  # ssdr <- summary(sdr)
}

# fit_mod(dat$tmb_data2[[4]], dat$tmb_pars2[[4]])

dat2 <- dat %>% 
  mutate(sdr = map2(tmb_data, tmb_pars, fit_mod),
         ssdr = map(sdr, summary),
         #also fit model to monthly effort data 
         sdr2 = map2(tmb_data2, tmb_pars2, fit_mod),
         ssdr2 = map(sdr2, summary))
saveRDS(dat2, here::here("generated_data", "model_fits", 
                         "combined_model_v3.RDS"))

map(dat2$sdr, print)
print(dat2$sdr[[4]])


## GENERATE PREDICTIONS --------------------------------------------------------

## Function to generate abundance predictions
gen_abund_pred <- function(catch, ssdr) {
  abund_pred <- ssdr[rownames(ssdr) %in% "agg_pred_abund", ] 
  # abund_b <- ssdr[rownames(ssdr) %in% "b1_j", ]
  # catch = FALSE because predictions are at region, not area scale
  pred_catch <- make_pred_dat(catch) %>% 
    mutate(reg_month = paste(region, month, sep = "_"))
  
  data.frame(pred_est = abund_pred[ , "Estimate"],
                        pred_se =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_catch %>% 
            select(region, month, reg_month) %>% 
            distinct(), .) %>% 
    mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
           pred_up = pred_est + (qnorm(0.975) * pred_se))
}

## Function to generate composition and stock-specific abundance predictions
gen_comp_pred <- function(infill_comp, ssdr) {
  comp_pred <- ssdr[rownames(ssdr) %in% "pred_probs", ]
  comp_abund_pred <- ssdr[rownames(ssdr) %in% "pred_abund_mg", ]
  
  stk_names <- unique(infill_comp$agg)
  pred_comp_dat <- make_pred_dat(infill_comp, catch = FALSE)
  n_preds <- nrow(pred_comp_dat)
  
  dum <- purrr::map_dfr(seq_along(stk_names), ~ pred_comp_dat)
  
  data.frame(stock = as.character(rep(stk_names, each = n_preds)),
                             pred_prob_est = comp_pred[ , "Estimate"],
                             pred_prob_se =  comp_pred[ , "Std. Error"]) %>% 
    cbind(dum, .) %>%
    mutate(pred_prob_low = pred_prob_est + (qnorm(0.025) * pred_prob_se),
           pred_prob_up = pred_prob_est + (qnorm(0.975) * pred_prob_se),
           comp_abund_est = comp_abund_pred[ , "Estimate"],
           comp_abund_se =  comp_abund_pred[ , "Std. Error"],
           comp_abund_low = comp_abund_est + (qnorm(0.025) * comp_abund_se),
           comp_abund_up = comp_abund_est + (qnorm(0.975) * comp_abund_se)) 
}



#predictions assuming mean effort
pred_dat <- dat2 %>% 
  mutate(abund_pred_ci = map2(catch, ssdr, gen_abund_pred),
         comp_pred_ci = map2(infill_comp, ssdr, gen_comp_pred))
# predictions assuming area-month specific effort (i.e. variable effort)
pred_dat2 <- dat2 %>% 
  mutate(abund_pred_ci = map2(catch, ssdr2, gen_abund_pred),
         comp_pred_ci = map2(infill_comp, ssdr2, gen_comp_pred))


## PLOT PREDICTIONS --------------------------------------------------------

file_path <- here::here("figs", "model_pred", "combined")

## Plot abundance
plot_abund <- function(dat, ylab) {
  ggplot() +
    geom_pointrange(data = dat, aes(x = as.factor(month), y = pred_est,
                                    ymin = pred_low, ymax = pred_up)) +
    labs(x = "", y = ylab) +
    ggsidekick::theme_sleek() +
    facet_wrap(~region)
  }

rec_catch <- plot_abund(pred_dat$abund_pred_ci[[1]], 
                        ylab = "Predicted\nMonthly Catch")
comm_catch <- plot_abund(pred_dat$abund_pred_ci[[2]], 
                         ylab = "Predicted\nDaily Catch")
rec_catch2 <- plot_abund(pred_dat2$abund_pred_ci[[1]], 
                        ylab = "Predicted\nMonthly Catch")
comm_catch2 <- plot_abund(pred_dat2$abund_pred_ci[[2]], 
                         ylab = "Predicted\nDaily Catch")
x.grob <- textGrob("Month", gp = gpar(col = "grey30", fontsize=10))

pdf(paste(file_path, "pred_abundance.pdf", sep = "/"))
cowplot::plot_grid(rec_catch, comm_catch, nrow = 2) %>% 
  arrangeGrob(., bottom = x.grob, 
              top = textGrob("Mean Effort", gp = gpar(fontsize=11))) %>% 
  grid.arrange
cowplot::plot_grid(rec_catch2, comm_catch2, nrow = 2) %>% 
  arrangeGrob(., bottom = x.grob, 
              top = textGrob("Variable Effort", gp = gpar(fontsize=11))) %>% 
  grid.arrange
dev.off()


## Composition predictions (effort doesn't impact predictions)
plot_comp <- function(comp_pred, raw_prop) {
  ggplot() +
    geom_point(data = raw_prop,
               aes(x = month, y = samp_g_ppn, fill = region),
               shape = 21, alpha = 0.4, position = position_dodge(0.6)) +
    geom_pointrange(data = comp_pred,
                    aes(x = month, y = pred_prob_est,
                        ymin = pred_prob_low, ymax = pred_prob_up,
                        fill = region),
                    shape = 21, size = 0.4, position = position_dodge(0.6)) +
    labs(y = "Probability", x = "Month") +
    facet_wrap(~ stock) +
    ggsidekick::theme_sleek()
}

pdf(paste(file_path, "composition.pdf", sep = "/"))
map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp)
dev.off()

# combined estimates of stock-specific CPUE
# note that effort is standardized differently between raw data and predictions
# so can't be easily compared
plot_ss_abund <- function(comp_pred) {
  ggplot() +
    # geom_point(data = pred_dat2$raw_abund[[1]], aes(x = month, y = cpue_g, fill = region),
    #            shape = 21, alpha = 0.3, position = position_dodge(0.6)) +
    geom_pointrange(data = comp_pred,
                    aes(x = month, y = comp_abund_est, ymin = comp_abund_low,
                        ymax = comp_abund_up, fill = region),
                    shape = 21, position = position_dodge(0.6)) +
    facet_wrap(~stock, ncol = 2, scales = "free_y") +
    scale_fill_viridis_d(option = "C", name = "Catch Region") +
    labs(x = "Month", y = "Predicted Catch") +
    ggsidekick::theme_sleek() +
    theme(legend.position="top")
}

pdf(paste(file_path, "stock-specific_abund_meanE.pdf", sep = "/"))
map(pred_dat$comp_pred_ci, plot_ss_abund)
dev.off()
pdf(paste(file_path, "stock-specific_abund_varyE.pdf", sep = "/"))
map(pred_dat2$comp_pred_ci, plot_ss_abund)
dev.off()


## DEFUNCT ##

## export plotting data for Rmd 
# list(catch = catch, pred_ci = pred_ci, raw_prop = raw_prop, raw_abund = raw_abund,
#      pred_catch = pred_catch) %>%
#   saveRDS(., here::here("generatedData", "model_fits", "pst_plot_list.RDS"))
# 
# 
# ## estimates of effort effects on catch
# n_betas <- length(coef(m1))
# eff_b <- abund_b[c(1, n_betas - 1, n_betas) , 1]
# 
# ggplot(catch) +
#   geom_point(aes(x = eff_z, y = catch), 
#              alpha = 0.2) +
#   # lims(x = c(0, 4)) +
#   stat_function(fun = function(x) exp(eff_b[1] + eff_b[2]*x + eff_b[3]*x^2)) 

