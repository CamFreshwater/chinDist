## Dirichlet model fit - observers only
# Nov. 19, 2020
# Fit dirichlet model to subset of stock composition data collected only by
# observers in rec fishery; used to evaluate impacts of using citizen scientists
# Otherwise duplicates fit_dirichlet_splines.R
# Most recent update: Dec. 9 2020

library(tidyverse)
library(TMB)
library(ggplot2)
library(mgcv)

# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(legal == "legal",
         !region %in% c("Queen Charlotte Sound")) %>%
  mutate(
    temp_strata = paste(month, region, sep = "_"),
    sample_id = paste(temp_strata, jDay, year, sep = "_"),
    min_m = case_when(
      region %in% c("S. Strait of Georgia", "Juan de Fuca Strait") ~ 4,
      region %in% c("N. Strait of Georgia", 
                    "Queen Charlotte and\nJohnstone Straits") ~ 6
    ),
    max_m = case_when(
      region %in% c("Juan de Fuca Strait", 
                    "N. Strait of Georgia",
                    "S. Strait of Georgia") ~ 9,
      region == "Queen Charlotte and\nJohnstone Straits" ~ 8
    )
  ) %>% 
  group_by(region) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>%
  ungroup() %>% 
  droplevels()

# subset fit including only data collected by creel observers
rec_obs <- rec %>% 
  filter(sampler == "Observer")

# helper function to calculate aggregate probs
calc_agg_prob <- function(grouped_data, full_data) {
  grouped_data %>% 
    summarize(agg_prob = sum(adj_prob), .groups = 'drop') %>% 
    arrange(sample_id, desc(agg_prob)) %>%
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
  if (raw_data$gear[1] == "sport") {
    temp %>%
      mutate(region = fct_relevel(region, "QCaJS", "JdFS", "NSoG", "SSoG"))
  } else {
    temp
  }
}

# combined tibble 
comp <- tibble(
  sample = c("gsi", "gsi", "gsi_obs", "gsi_obs"),
  fishery = rep("sport", times = 4),
  grouping = rep(c("pst", "can"), 2),
  dataset = paste(sample, fishery, grouping, sep = "_"),
  # incorporate raw_comp data
  raw_data = list(rec, rec, rec_obs, rec_obs)
) %>% 
  mutate(
    grouping_col = case_when(
      grouping == "pst" ~ "pst_agg",
      grouping == "can" ~ "reg1"
    ),
    comp_long = pmap(list(grouping_col, raw_data), .f = clean_comp)
  ) 


## PREP INPUTS ----------------------------------------------------------------

# add tmb data assuming average effort predictions
tmb_list <- purrr::map(comp$comp_long, 
                       .f = stockseasonr::gen_tmb,
                       random_walk = TRUE,
                       model_type = "composition")

# join all data into a tibble then generate model inputs
dat <- comp %>% 
  mutate(
    comp_wide = purrr::map(tmb_list, "comp_wide"),
    tmb_data = purrr::map(tmb_list, "data"),
    tmb_pars = purrr::map(tmb_list, "pars"),
    tmb_map = purrr::map(tmb_list, "tmb_map"),
    pred_dat_comp = purrr::map(tmb_list, ~ .$pred_dat_comp)
  ) 


# FIT --------------------------------------------------------------------------

# optim_fix_inits_vec = c(FALSE, FALSE)
# 
# dat2 <- dat %>%
#   mutate(sdr = purrr::pmap(list(comp_dat = .$comp_long,
#                                 optim_fix_inits = FALSE),
#                            .f = stockseasonr::fit_stockseason,
#                            random_walk = TRUE,
#                            model_type = "composition",
#                            silent = FALSE))
# dat2$ssdr <- purrr::map(dat2$sdr, summary)
# 
# saveRDS(dat2 %>% select(-sdr),
#         here::here("generated_data", "model_fits", 
#                    "composition_model_dir_supp.RDS"))

dat2 <- readRDS(here::here("generated_data", "model_fits", 
                           "composition_model_dir_supp.RDS"))


## GENERATE OBS AND PREDICTIONS ------------------------------------------------

source(here::here("R", "functions", "plot_cleaning_functions.R"))
 
comp_pred_dat <- dat2 %>%
  mutate(
    comp_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr),
                                         .f = gen_comp_pred,
                                         comp_only = TRUE)),
    raw_prop = purrr::map2(comp_long, comp_wide, make_raw_prop_dat)
  ) %>%
  select(dataset:comp_wide, comp_pred_ci, raw_prop) %>%
  mutate(
    comp_pred_ci = map2(grouping_col, comp_pred_ci, .f = stock_reorder),
    raw_prop = map2(grouping_col, raw_prop, .f = stock_reorder)
  )


## PLOT PREDICTIONS ------------------------------------------------------------

source(here::here("R", "functions", "plot_predictions_splines.R"))

#color palette
pal <- readRDS(here::here("data", "disco_color_pal.RDS"))

#pfma map (used to steal legend)
pfma_map <- readRDS(here::here("generated_data", "pfma_map.rds"))

# generic settings
file_path <- here::here("figs", "model_pred", "dirichlet_only")

# Composition spline prediction
comp_plots <- map2(comp_pred_dat$comp_pred_ci, comp_pred_dat$raw_prop, plot_comp,
                   raw = TRUE)
comp_plots_fix <- map2(comp_pred_dat$comp_pred_ci, comp_pred_dat$raw_prop, plot_comp,
                       raw = TRUE, facet_scales = "fixed")

pdf(paste(file_path, "composition_splines_supp.pdf", sep = "/"))
comp_plots
dev.off()


## Composition stacked ribbon prediction
combine_data <- function(in_col = c("pst_agg", "reg1")) {
  dum <- comp_pred_dat %>%
    filter(grouping_col == in_col) %>%
    select(dataset, grouping_col, comp_pred_ci) %>%
    unnest(., cols = c(comp_pred_ci)) %>% 
    mutate(region_c = fct_relevel(region_c, 
                                  "Queen Charlotte and\nJohnstone Straits", 
                                  "Juan de Fuca Strait", 
                                  "N. Strait of Georgia", 
                                  "S. Strait of Georgia"),
           dataset_label = case_when(
             grepl("obs", dataset) ~ "Observers Only",
             TRUE ~ "Observers\nand Volunteers")
           ) 
  
  # consolidate Col Spring to reduce total number of categories and improve 
  # readability
  if (in_col == "pst_agg") {
    dum <- dum %>% 
      mutate(
        stock = fct_recode(stock, "CR-spring" = "CR-lower_sp", 
                           "CR-spring" = "CR-upper_sp")) %>% 
      group_by(stock, region, region_c, month_n, dataset_label) %>% 
      summarize(pred_prob_est = sum(pred_prob_est), .groups = "drop")
  }
  
  return(dum)
}

pst_comp_pred <- combine_data(in_col = "pst_agg") 
pst_area <- plot_comp_stacked(pst_comp_pred, grouping_col = "pst_agg", 
                              palette_name = "sunset") +
  facet_grid(region_c ~ dataset_label)
can_comp_pred <- combine_data(in_col = "reg1") 
can_area <- plot_comp_stacked(can_comp_pred, grouping_col = "reg1", 
                              palette_name = "midnight") +
  facet_grid(region_c ~ dataset_label)

pdf(paste(file_path, "composition_stacked_supp.pdf", sep = "/"))
pst_area
can_area
dev.off()

#save plot as list to export to fit_dirichlet_splines for printing
comp_stacked_supp_list <- list(pst_area, can_area)
saveRDS(comp_stacked_supp_list, 
        here::here("generated_data", "comp_stacked_supp_list.RDS"))
