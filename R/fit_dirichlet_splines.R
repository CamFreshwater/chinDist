## Dirichlet model fit
# Oct 29, 2020
# Fit dirichlet model to stock composition data
# Duplicates composition component of combined model (fit_combined_splines.R)
# Ran independently to increase monthly span

library(tidyverse)
library(ggplot2)

# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(legal == "legal",
         !region %in% c("Queen Charlotte Sound")) %>%
  mutate(
    temp_strata = paste(month, region, sep = "_"),
    sample_id = paste(temp_strata, jDay, year, sep = "_"),
    min_m = case_when(
      region %in% c("N. Strait of Georgia", "S. Strait of Georgia") ~ 1,
      region == "Queen Charlotte and\nJohnstone Straits" ~ 6,
      region == "Juan de Fuca Strait" ~ 3
    ),
    max_m = case_when(
      region %in% c("Juan de Fuca Strait", #"S. Strait of Georgia", 
                    "N. Strait of Georgia") ~ 9,
      region == "S. Strait of Georgia" ~ 12,
      region == "Queen Charlotte and\nJohnstone Straits" ~ 8
    )
  ) %>% 
  group_by(region) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>%
  ungroup() %>% 
  droplevels()

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
    comp_long = pmap(list(grouping_col, raw_data), .f = clean_comp)
  ) 


## PREP INPUTS ----------------------------------------------------------------

# add tmb data assuming average effort predictions
# note that these outputs are generated automatically in fit_stockseason
# but aren't retained and are necessary for the current suite of plotting 
# functions
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
# 
# #use fix optim inits except rec_pst (convergence issues)
optim_fix_inits_vec = c(TRUE, FALSE, TRUE, TRUE)
# optim_fix_inits_vec = c(TRUE, TRUE)

dat2 <- dat %>%
  # filter(!fishery == "sport") %>%
  mutate(sdr = purrr::pmap(list(comp_dat = .$comp_long,
                                optim_fix_inits = optim_fix_inits_vec),
                           .f = stockseasonr::fit_stockseason,
                           random_walk = TRUE,
                           model_type = "composition",
                           silent = FALSE))
dat2$ssdr <- purrr::map(dat2$sdr, summary)

# dat2 <- rbind(dat4[1,], dat3[1,], dat4[2,], dat3[2, ])

saveRDS(dat2 %>% select(-sdr),
        here::here("generated_data", "model_fits", "composition_model_dir.RDS"))
dat2 <- readRDS(here::here("generated_data", "model_fits", "composition_model_dir.RDS"))

## GENERATE OBS AND PREDICTIONS ------------------------------------------------

# source(here::here("R", "functions", "plot_cleaning_functions.R"))
# 
# comp_pred_dat <- dat2 %>%
#   mutate(
#     comp_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr),
#                                          .f = gen_comp_pred,
#                                          comp_only = TRUE)),
#     raw_prop = purrr::map2(comp_long, comp_wide, make_raw_prop_dat)
#   ) %>%
#   select(dataset:comp_wide, comp_pred_ci, raw_prop) %>%
#   mutate(
#     comp_pred_ci = map2(grouping_col, comp_pred_ci, .f = stock_reorder),
#     raw_prop = map2(grouping_col, raw_prop, .f = stock_reorder)
#   )
# saveRDS(comp_pred_dat, here::here("generated_data",
#                              "composition_model_predictions.RDS"))
comp_pred_dat <- readRDS(here::here("generated_data",
                               "composition_model_predictions.RDS"))


## PLOT PREDICTIONS ------------------------------------------------------------

source(here::here("R", "functions", "plot_predictions_splines.R"))

#color palette
pal <- readRDS(here::here("generated_data", "disco_color_pal.RDS"))

#pfma map (used to steal legend)
pfma_map <- readRDS(here::here("generated_data", "pfma_map.rds"))

# generic settings
file_path <- here::here("figs", "model_pred", "dirichlet_only")

# Composition spline prediction
comp_plots <- map2(comp_pred_dat$comp_pred_ci, comp_pred_dat$raw_prop, plot_comp,
                   raw = TRUE)
comp_plots_fix <- map2(comp_pred_dat$comp_pred_ci, comp_pred_dat$raw_prop, plot_comp,
                       raw = TRUE, facet_scales = "fixed")

pdf(paste(file_path, "composition_splines.pdf", sep = "/"))
comp_plots
comp_plots_fix
dev.off()


## Composition stacked ribbon prediction
combine_data <- function(in_col = c("pst_agg", "reg1")) {
  comp_pred_dat %>%
    filter(grouping_col == in_col) %>%
    select(dataset, grouping_col, comp_pred_ci) %>%
    unnest(., cols = c(comp_pred_ci)) %>% 
    mutate(region_c = fct_relevel(region_c, "NWVI", "SWVI", 
                                  "Queen Charlotte and\nJohnstone Straits", 
                                  "Juan de Fuca Strait", 
                                  "N. Strait of Georgia", "S. Strait of Georgia"))
}
pst_comp_pred <- combine_data(in_col = "pst_agg") 
pst_area <- plot_comp_stacked(pst_comp_pred, grouping_col = "pst_agg")
can_comp_pred <- combine_data(in_col = "reg1") 
can_area <- plot_comp_stacked(can_comp_pred, grouping_col = "reg1")

pdf(paste(file_path, "composition_stacked.pdf", sep = "/"))
pst_area
can_area
dev.off()


# MS figs
png(here::here("figs", "ms_figs", "comp_pst_stacked.png"), res = 400, units = "in",
    height = 6, width = 7.5)
pst_area
dev.off()

png(here::here("figs", "ms_figs", "comp_can_stacked.png"), res = 400, units = "in",
    height = 6, width = 7.5)
can_area
dev.off()

# Supplementary figs
png(here::here("figs", "ms_figs", "comp_ext_pst_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
comp_plots_fix[[1]]
dev.off()

png(here::here("figs", "ms_figs", "comp_ext_pst_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
comp_plots_fix[[2]]
dev.off()

png(here::here("figs", "ms_figs", "comp_ext_can_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7)
comp_plots_fix[[3]]
dev.off()

png(here::here("figs", "ms_figs", "comp_ext_can_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7)
comp_plots_fix[[4]]
dev.off()

# supplementary composition analysis from fit_dirichlet_splines_supp.R
comp_stacked_supp_list <- readRDS(here::here("generated_data", 
                                             "comp_stacked_supp_list.RDS"))
png(here::here("figs", "ms_figs", "comp_pst_stacked_supp.png"), res = 400, units = "in",
    height = 6.5, width = 7)
comp_stacked_supp_list[[1]]
dev.off()
png(here::here("figs", "ms_figs", "comp_can_stacked_supp.png"), res = 400, units = "in",
    height = 6.5, width = 7)
comp_stacked_supp_list[[2]]
dev.off()
