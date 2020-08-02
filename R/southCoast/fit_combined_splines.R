## Combined model fits
# July 20, 2020
# Fit combined multionomial/tweedie or multinomial/nb model to 
# stock composition and abundance data
# aggregate probability reflects summed probabilities of a given region of 
# origin for a given individual
# Modified from fit_combined.R to include splines

library(tidyverse)
library(TMB)
library(ggplot2)
library(grid)
library(gridExtra)
library(mgcv)
library(scales)

# CLEAN CATCH -----------------------------------------------------------------

# clean function 
clean_catch <- function(dat) {
  dat %>% 
    filter(!is.na(eff),
           !eff == "0") %>% 
    mutate(region_c = as.character(region),
           region = factor(abbreviate(region, minlength = 4)),
           area_n = as.numeric(area),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff))
           ) %>% 
    arrange(region, month) 
}

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  # calculate monthly catch and effort 
  group_by(catchReg, area, month, year) %>% 
  summarize(catch = sum(catch), 
            eff = sum(boatDays)) %>% 
  ungroup() %>% 
  rename(region = catchReg) %>% 
  clean_catch(.) %>% 
  # drop inside areas where seasonal catches not available
  filter(!area_n < 100) %>% 
  droplevels() %>% 
  select(region, region_c, area, area_n, month, month_n, year, catch, eff, 
         eff_z)
  
#recreational catch data
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "monthlyCatch_rec.RDS")) %>% 
  #drop sublegal fish and regions without genetics
  filter(legal == "legal",
         !region %in% c("NWVI", "SWVI")) %>% 
  #group by subarea to get rid of adipose and released legal duplicates
  group_by(month, month_n, year, area, subarea, region, legal) %>% 
  mutate(subarea_catch = sum(mu_catch),
         subarea_eff = mean(mu_boat_trips)) %>% 
  #group by area to be consistent with commercial data
  group_by(month, month_n, year, area, region, legal) %>% 
  summarize(catch = sum(subarea_catch),
            eff = sum(subarea_eff)) %>% 
  ungroup() %>% 
  clean_catch(.) %>% 
  # drop months with minimal catch estimates
  filter(!month_n < 3,
         !month_n > 10) %>% 
  # drop areas with fewer than 10 datapoints
  group_by(area_n) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(!n < 10) %>% 
  droplevels() %>% 
  select(region, region_c, area, area_n, month, month_n, year, catch, eff, 
         eff_z)

# export areas to make maps
# areas_retained <- c(unique(rec_catch$area_n), unique(comm_catch$area_n))
# saveRDS(areas_retained, 
#         here::here("data", "gsiCatchData", "pfma", "areas_to_plot.RDS"))


# CLEAN GENETICS  --------------------------------------------------------------

# recreational data
rec <-  readRDS(here::here("data", "gsiCatchData", "rec", 
                           "recIndProbsLong.rds")) %>% 
  filter(legal == "legal",
         !month_n < min(rec_catch$month_n),
         !month_n > max(rec_catch$month_n)) %>%
  mutate(sample_id = paste(temp_strata, jDay, year, sep = "_"))

# check which WCVI stocks are in Johnstone Strait
# rec %>% 
#   filter(region == "Johnstone Strait") %>% 
#   group_by(sample_id, stock, Region1Name, area) %>% 
#   summarize(lump_prob = sum(adj_prob)) %>% 
#   filter(Region1Name == "WCVI",
#          !lump_prob == 0) %>%
#   group_by(stock) %>% 
#   filter(!lump_prob < 4) %>% 
#   ungroup() %>% 
#   ggplot(.) +
#   geom_boxplot(aes(x = stock, y = lump_prob)) +
#   facet_wrap(~area)

# GSI samples in commercial database associated with Taaq fishery
taaq <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                            "taaq_summary.csv"), stringsAsFactors = F) %>% 
  filter(drop == "y") %>% 
  mutate(temp_strata2 = paste(month, region, year, sep = "_"))
# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds")) %>% 
  # drop month-region-years where catch data are missing or where samples came 
  # from Taaq fishery, and no comm fishery was active in the same region
  mutate(sample_id = paste(temp_strata, jDay, year, sep = "_"),
         region_c = as.character(region), 
         temp_strata2 = paste(month, region, year, sep = "_")) %>% 
  filter(!temp_strata2 %in% taaq$temp_strata2) %>%
  select(-temp_strata2) %>% 
  droplevels()

# helper function to calculate aggregate probs
calc_agg_prob <- function(grouped_data, full_data) {
  grouped_data %>% 
    summarize(agg_prob = sum(adj_prob)) %>% 
    arrange(sample_id, desc(agg_prob)) %>%
    ungroup() %>% 
    distinct() %>% 
    left_join(.,
              full_data %>% 
                #important to subset appropriately to remove individual traits
                #that lead to duplicates
                select(sample_id, region, year, month, gear, month_n) %>% 
                distinct(),
              by = "sample_id")
}

# helper function to pool non-focal stocks
pool_aggs <- function(full_data) {
  full_data %>% 
    mutate(
      # pst_agg = case_when(
      #   grepl("CR-", pst_agg) ~ "CR",
      #   grepl("CST", pst_agg) ~ "CA/OR/WA",
      #   pst_agg %in% c("NBC_SEAK", "WCVI") ~ "BC-coast",
      #   pst_agg %in% c("PSD", "SOG") ~ "SalSea",
      #   TRUE ~ pst_agg
      # ),
      reg1 = case_when(
        Region1Name %in% c("Fraser_Fall", "ECVI", "Fraser_Summer_4.1", 
                           "Fraser_Spring_4.2", "Fraser_Spring_5.2", 
                           "Fraser_Summer_5.2", "WCVI", "SOMN") ~ Region1Name,
        TRUE ~ "Other"
      )
    )
}

# helper function to use above for gsi and cwt samples (cwt ignored for now)
clean_comp <- function(grouping_col, raw_data, ...) {
  temp <- raw_data %>% 
    pool_aggs() %>% 
    rename(agg = grouping_col) %>%
    group_by(sample_id, agg) %>%
    calc_agg_prob(., raw_data) %>% 
    group_by(region, month, year) %>%
    mutate(nn = sum(agg_prob)) %>%
    #remove strata with less than 10 individuals total
    filter(!nn < 10) %>%
    ungroup() %>%
    droplevels() %>%
    mutate(region_c = as.character(region),
           region = factor(abbreviate(region, minlength = 4))) %>% 
    select(sample_id, gear, region, region_c, year, month, month_n, agg,
           agg_prob, nn) %>%
    distinct()
}

# temp %>% select(region, year, month, nn) %>% distinct() %>% pull(nn) %>% sum()
# temp %>% select(region, year, month) %>% distinct() %>% nrow()

## combine genetics in tibble 
full_dat <- tibble(
  sample = rep("gsi", 4),
  fishery = rep(c("troll", "sport"), times = 2),
  grouping = rep(c(rep("pst", 2), rep("can", 2)), 1),
  dataset = paste(sample, fishery, grouping, sep = "_"),
  # incorporate raw_comp data
  raw_data = list(comm, rec, comm, rec)
) %>% 
  mutate(
    grouping_col = case_when(
      grouping == "pst" ~ "pst_agg",
      grouping == "can" ~ "reg1"
    ),
    comp_long = pmap(list(grouping_col, raw_data), .f = clean_comp),
    catch_data = list(comm_catch, rec_catch, comm_catch, rec_catch)
  )

# Check coverage of recreational fishery
tt <- full_dat %>% 
  filter(dataset == "gsi_sport_pst")
# gsi
tt$comp_long[[1]] %>% 
  select(year, month, region, nn) %>%
  distinct() %>% 
  group_by(year, month, region) %>%
  mutate(strata_n = sum(nn)) %>% 
  arrange(region, month, year) %>% 
  print(n = Inf)

# catch 
tt$catch_data[[1]] %>%
  group_by(month, region) %>% 
  tally() %>% 
  arrange(region, month) %>% 
  print(n = Inf)

# How many samples total?
f <- function(x) {
  x %>% 
    select(month, year, region, nn) %>% 
    distinct() %>% 
    pull(nn) %>% 
    sum()
}
f(full_dat$comp_long[[1]])
f(full_dat$comp_long[[2]])

## PREP MODEL INPUTS -----------------------------------------------------------

catch_dat <- full_dat$catch_data[[2]]
comp_dat <- full_dat$comp_long[[2]]

# function to generate TMB data inputs from wide dataframes
gen_tmb <- function(catch_dat, comp_dat) {
  ## catch data
  yr_vec_catch <- as.numeric(as.factor(as.character(catch_dat$year))) - 1
  
  #generate model matrix based on GAM with spline type a function of months 
  #retained
  months1 <- unique(catch_dat$month_n)
  spline_type <- ifelse(max(months1) == 12, "cc", "tp")
  n_knots <- ifelse(max(months1) == 12, 4, 3)
  m1 <- gam(catch ~ 
              area + s(month_n, bs = spline_type, k = n_knots, by = area) +
              s(eff_z, bs = "tp", k = 4),
            # knots = list(month_n = c(min(months1), max(months1))),
            data = catch_dat,
            family = nb)
  fix_mm_catch <- predict(m1, type = "lpmatrix")

  # generic data frame that is skeleton for both comp and catch predictions
  if (comp_dat$gear[1] == "sport") {
    pred_dat <- group_split(comp_dat, region) %>% 
      map_dfr(., function (x) {
        expand.grid(
          month_n = seq(min(x$month_n), 
                        max(x$month_n),
                        by = 0.1),
          region = unique(x$region)
        )
      }) 
  } else {
    pred_dat <- expand.grid(
      month_n = seq(min(comp_dat$month_n), 
                    max(comp_dat$month_n),
                    by = 0.1),
      region = unique(comp_dat$region)
    )
  }
  # list of strata from composition dataset to retain in catch predictive model 
  # to ensure aggregate abundance can be calculated
  comp_strata <- unique(paste(pred_dat$month_n, pred_dat$region, sep = "_"))
  
  # make predictive model matrix including null values for effort
  pred_dat_catch <- expand.grid(
    month_n = seq(min(catch_dat$month_n),
                  max(catch_dat$month_n),
                  by = 0.1),
    area = unique(catch_dat$area),
    eff_z = 0
  ) %>%
    left_join(., 
              catch_dat %>% select(region, area) %>% distinct(), 
              by = "area") %>% 
    mutate(strata = paste(month_n, region, sep = "_")) %>% 
    filter(strata %in% comp_strata) %>% 
    #convoluted fuckery to ensure ordering is correct for key
    mutate(month = as.factor(month_n),
           order = row_number(),
           reg_month = paste(region, month, sep = "_"),
           reg_month_f = fct_reorder(factor(reg_month), order)) %>% 
    select(-order, -reg_month, -strata) %>%
    distinct()
  pred_mm_catch <- predict(m1, pred_dat_catch, type = "lpmatrix")
  
  # construct factor key for regional aggregates associated with areas
  grouping_vec <- as.numeric(pred_dat_catch$reg_month_f) - 1
  grouping_key <- unique(grouping_vec)

  
  ## composition data
  gsi_wide <- comp_dat %>% 
    pivot_wider(., names_from = agg, values_from = agg_prob) %>%
    mutate_if(is.numeric, ~replace_na(., 0.000001))
  
  obs_comp <- gsi_wide %>% 
    select(-c(sample_id:nn)) %>% 
    as.matrix() 
  
  yr_vec_comp <- as.numeric(gsi_wide$year) - 1
  
  months2 <- unique(gsi_wide$month_n)
  spline_type <- ifelse(max(months2) == 12, "cc", "tp")
  n_knots <- ifelse(max(months2) == 12, 4, 3)
  # response variable doesn't matter, since not fit
  m2 <- gam(rep(0, length.out = nrow(gsi_wide)) ~ 
              region + s(month_n, bs = spline_type, k = n_knots, by = region), 
            # knots = list(month_n = c(min(months2), max(months2))),
            data = gsi_wide)
  fix_mm_comp <- predict(m2, type = "lpmatrix")
  
  pred_mm_comp <- predict(m2, pred_dat, type = "lpmatrix")

  #input data
  data = list(
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
  
  #input parameter initial values
  pars = list(
    #abundance parameters
    b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
    log_phi = log(1.5),
    z1_k = rep(0, length(unique(yr_vec_catch))),
    log_sigma_zk1 = log(0.25),
    #composition parameters
    b2_jg = matrix(0, 
                   nrow = ncol(fix_mm_comp), 
                   ncol = ncol(obs_comp)),
    z2_k = rep(0, times = length(unique(yr_vec_comp))),
    log_sigma_zk2 = log(0.75)
  )
  list("data" = data, "pars" = pars, "comp_wide" = gsi_wide,
       "pred_dat_catch" = pred_dat_catch, "pred_dat_comp" = pred_dat)
}

#add tmb data assuming average effort predictions (for month specific look at
# non-spline code)
tmb_list <- map2(full_dat$catch_data, full_dat$comp_long, gen_tmb)


# join all data into a tibble then generate model inputs
dat <- full_dat %>% 
  mutate(
    comp_wide = purrr::map(tmb_list, "comp_wide"),
    tmb_data = purrr::map(tmb_list, "data"),
    tmb_pars = purrr::map(tmb_list, "pars"),
    pred_dat_catch = purrr::map(tmb_list, ~ .$pred_dat_catch),
    pred_dat_comp = purrr::map(tmb_list, ~ .$pred_dat_comp)
  )


## FIT MODELS ------------------------------------------------------------------

compile(here::here("src", "nb_dirichlet_1re.cpp"))
dyn.load(dynlib(here::here("src", "nb_dirichlet_1re")))

fit_mod <- function(tmb_data, tmb_pars, nlminb_loops = 2) {
  obj <- MakeADFun(tmb_data, tmb_pars, 
                   # dat$tmb_data[[4]], dat$tmb_pars[[4]],
                   random = c("z1_k", "z2_k"), 
                   DLL = "nb_dirichlet_1re")
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  if (nlminb_loops > 1) {
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr)
    }
  }
  sdreport(obj)
  # sdr <- sdreport(obj)
  # ssdr <- summary(sdr)
}

# join all data into a tibble then generate model inputs
dat2 <- dat %>%
  mutate(sdr = purrr::map2(tmb_data, tmb_pars, fit_mod, nlminb_loops = 2),
         ssdr = purrr::map(sdr, summary)
         )
saveRDS(dat2 %>% select(-sdr), 
        here::here("generated_data", "model_fits", "combined_model_dir.RDS"))

dat2 <- readRDS(here::here("generated_data", "model_fits",
                           "combined_model_dir.RDS"))

## GENERATE OBS AND PREDICTIONS ------------------------------------------------

# function to calculate raw proportions based on observed composition data
make_raw_prop_dat <- function(comp_long, comp_wide) {
  n_stks <- length(unique(comp_long$agg))
  comp_wide %>%
    pivot_longer((ncol(.) - n_stks + 1):ncol(.), 
                 names_to = "agg", values_to = "agg_prob") %>% 
    group_by(region, region_c, month_n, year, agg) %>%
    summarize(samp_g = sum(agg_prob)) %>%
    group_by(region, region_c, month_n, year) %>%
    mutate(samp_total = sum(samp_g)) %>% 
    ungroup() %>% 
    mutate(samp_g_ppn = samp_g / samp_total,
           stock = fct_reorder(agg, desc(samp_g_ppn))) 
}
  
# function to calculate raw stock-specific abundance observations
# data_type necessary for monthly vs. daily predictions
make_raw_abund_dat <- function(catch, raw_prop, fishery) {
  if (fishery == "sport") {
    dum <-  catch %>% 
      group_by(region, region_c, month_n, year) %>%
      summarize(
        #sum rec because subareas each have one monthly value
        #mean commercial because daily values for each area
        cum_catch = as.numeric(sum(catch)),
        cum_effort = as.numeric(sum(eff))) %>% 
      ungroup()
  }
  if (fishery == "troll") {
    dum <- catch %>% 
      group_by(region, region_c, area, month_n, year) %>% 
      # calculate daily catch within month first
      summarize(
        mu_catch = mean(catch),
        mu_effort = mean(eff)
      ) %>% 
      group_by(region, region_c, month_n, year) %>% 
      summarize(
        cum_catch = as.numeric(sum(mu_catch)),
        cum_effort = as.numeric(sum(mu_effort))
      ) %>% 
      ungroup()
  }
  dum  %>%
    mutate(agg_cpue = cum_catch / cum_effort,
           month = as.character(month_n)) %>%
    left_join(., raw_prop, by = c("region", "month_n", "year")) %>%
    mutate(catch_g = samp_g_ppn * cum_catch,
           cpue_g = samp_g_ppn * agg_cpue,
           region = as.factor(region)) %>%
    filter(!is.na(stock))
}

## Function to generate aggregate abundance predictions (use comp predictions
# because aggregate)
gen_abund_pred <- function(comp_long, pred_dat_comp, ssdr) {
  abund_pred <- ssdr[rownames(ssdr) %in% "agg_pred_abund", ] 
  pred_dat <- pred_dat_comp %>% 
    left_join(., comp_long %>% select(region, region_c) %>% distinct(), 
              by = "region")
  
  data.frame(pred_est = abund_pred[ , "Estimate"],
                        pred_se =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_dat, .) %>% 
    mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
           pred_up = pred_est + (qnorm(0.975) * pred_se))
}

## Function to generate composition and stock-specific abundance predictions
gen_comp_pred <- function(comp_long, pred_dat_comp, ssdr) {
  comp_pred <- ssdr[rownames(ssdr) %in% "pred_pi_prop", ]
  comp_abund_pred <- ssdr[rownames(ssdr) %in% "pred_abund_mg", ]
  
  stk_names <- unique(comp_long$agg)
  n_preds <- nrow(pred_dat_comp)
  
  pred_dat <- pred_dat_comp %>% 
    left_join(., comp_long %>% select(region, region_c) %>% distinct(), 
              by = "region") 
  dum <- purrr::map_dfr(seq_along(stk_names), ~ pred_dat)
  
  data.frame(stock = as.character(rep(stk_names, each = n_preds)),
             pred_prob_est = comp_pred[ , "Estimate"],
             pred_prob_se =  comp_pred[ , "Std. Error"]) %>% 
    cbind(dum, .) %>%
    mutate(pred_prob_low = pmax(pred_prob_est + (qnorm(0.025) * pred_prob_se),
                                0),
           pred_prob_up = pmin(pred_prob_est + (qnorm(0.975) * pred_prob_se), 
                               1),
           comp_abund_est = as.vector(comp_abund_pred[ , "Estimate"]),
           comp_abund_se =  as.vector(comp_abund_pred[ , "Std. Error"]),
           comp_abund_low = pmax(comp_abund_est + (qnorm(0.025) * comp_abund_se),
                                 0),
           comp_abund_up = comp_abund_est + (qnorm(0.975) * comp_abund_se)) 
}

## Function to reorder stocks for plotting
stock_reorder <- function(key, comp_data) {
  if (key == "pst_agg") {
    out <- comp_data %>% 
      mutate(stock =  fct_relevel(stock, "CA_ORCST", "CR-sp&su", "CR-bright", 
                                  "CR-tule", "WACST", "PSD", "SOG", "FR-early", 
                                  "FR-late", "WCVI", "NBC_SEAK"),
             stock = fct_recode(stock, "WA-coast" = "WACST", 
                                "CA/OR-coast" = "CA_ORCST", 
                                "NBC/SEAK" = "NBC_SEAK"))
  }
  if (key == "reg1") {
    out <- comp_data %>% 
      mutate(stock =  fct_relevel(stock, "Fraser_Spring_4.2", 
                                  "Fraser_Spring_5.2", "Fraser_Summer_4.1", 
                                  "Fraser_Summer_5.2", "Fraser_Fall", "SOMN", 
                                  "ECVI", "WCVI", "Other"
      ))
  }
  return(out)
}


pred_dat <- dat2 %>% 
  mutate(
    #predictions assuming mean effort
    abund_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr), 
                                          .f = gen_abund_pred)),
    comp_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr), 
                        .f = gen_comp_pred)),
    raw_prop = purrr::map2(comp_long, comp_wide, make_raw_prop_dat),
    raw_abund = pmap(list(catch_data, raw_prop, fishery), 
                     .f = make_raw_abund_dat)) %>% 
  select(dataset:catch_data, abund_pred_ci:raw_abund) %>% 
  mutate(
    comp_pred_ci = map2(grouping_col, comp_pred_ci, .f = stock_reorder),
    raw_prop = map2(grouping_col, raw_prop, .f = stock_reorder)
  )
saveRDS(pred_dat, here::here("generated_data", 
                             "combined_model_dir_predictions.RDS"))
# pred_dat <- readRDS(here::here("generated_data",
#                                "combined_model_predictions.RDS"))


## PLOT PREDICTIONS ------------------------------------------------------------

# source file for plotting functions
source(here::here("R", "functions", "plot_predictions_splines.R"))

#color palette
pal <- readRDS(here::here("generated_data", "color_pal.RDS"))

#pfma map (used to steal legend)
pfma_map <- readRDS(here::here("generated_data", "pfma_map.rds"))

# generic settings
file_path <- here::here("figs", "model_pred", "combined")

## Plot abundance
comm1 <- pred_dat %>% 
  filter(dataset == "gsi_troll_pst") 
comm_abund <- plot_abund(comm1$abund_pred_ci[[1]], 
                         ylab = "Predicted Daily\nCatch Index") +
  annotate("text", x = -Inf, y = Inf, label = "a)", hjust = -1, vjust = 2)
rec1 <- pred_dat %>% 
  filter(dataset == "gsi_sport_pst") 
rec_abund <- plot_abund(rec1$abund_pred_ci[[1]], 
                        ylab = "Predicted Monthly\nCatch Index") +
  annotate("text", x = -Inf, y = Inf, label = "b)", hjust = -1, vjust = 2)

combo_abund <- cowplot::plot_grid(comm_abund, rec_abund, nrow = 2) %>% 
  arrangeGrob(., 
              bottom = textGrob("Month", 
                                gp = gpar(col = "grey30", fontsize=10))) %>% 
  grid.arrange()
# add legend from PFMA map
combo_abund2 <- cowplot::plot_grid(
  cowplot::get_legend(pfma_map),
  combo_abund,
  ncol=1, rel_heights=c(.1, .9)
)

pdf(paste(file_path, "pred_abundance_splines.pdf", sep = "/"))
combo_abund2
dev.off()

## Composition prediction
comp_plots <- map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp, 
                   raw = TRUE)
comp_plots_fix <- map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp, 
                       raw = TRUE, facet_scales = "fixed")
comp_plots_fix2 <- map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp, 
                        raw = FALSE, facet_scales = "fixed")

pdf(paste(file_path, "composition_splines.pdf", sep = "/"))
comp_plots
comp_plots_fix
dev.off()

## Stock-specific CPUE predictions
ss_abund_plots_free <- map2(pred_dat$comp_pred_ci, pred_dat$raw_abund, 
                            plot_ss_abund, raw = FALSE)
ss_abund_plots_fix <- map2(pred_dat$comp_pred_ci, pred_dat$raw_abund, 
                           plot_ss_abund,  raw = FALSE, facet_scales = "fixed")
pdf(paste(file_path, "stock-specific_abund_splines.pdf", sep = "/"))
ss_abund_plots_free
ss_abund_plots_fix
dev.off()


## Manuscript figures   
png(here::here("figs", "ms_figs", "abund_pred.png"), res = 400, units = "in",
    height = 5, width = 5)
combo_abund2
dev.off()

png(here::here("figs", "ms_figs", "comp_pst_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
comp_plots_fix[[1]]
dev.off()

png(here::here("figs", "ms_figs", "comp_pst_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
comp_plots_fix[[2]]
dev.off()

png(here::here("figs", "ms_figs", "comp_can_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7)
comp_plots_fix[[3]]
dev.off()

png(here::here("figs", "ms_figs", "comp_can_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7)
comp_plots_fix[[4]]
dev.off()

png(here::here("figs", "ms_figs", "abund_pst_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
ss_abund_plots_free[[1]] +
  labs(y = "Predicted Daily Catch Index")
dev.off()

png(here::here("figs", "ms_figs", "abund_pst_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
ss_abund_plots_free[[2]] +
  labs(y = "Predicted Monthly Catch Index")
dev.off()

png(here::here("figs", "ms_figs", "abund_can_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7)
ss_abund_plots_free[[3]] +
  labs(y = "Predicted Daily Catch Index")
dev.off()

png(here::here("figs", "ms_figs", "abund_can_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7)
ss_abund_plots_free[[4]] +
  labs(y = "Predicted Monthly Catch Index")
dev.off()