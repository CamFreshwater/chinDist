## Combined model fits
# Sep 25, 2020
# Fit combined dirichlet/nb model to stock composition and abundance data
# aggregate probability reflects summed probabilities of a given region of 
# origin for a given individual
# Constrain to relatively few months due to limited availability of rec C/E 
# data
# Modified from fit_combined.R to include splines

library(tidyverse)
library(TMB)
library(ggplot2)
library(grid)
library(gridExtra)
library(mgcv)
library(scales)

# IMPORT CATCH -----------------------------------------------------------------

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "month_area_commCatch.rds")) 

#recreational catch data - sampling unit is area-month-year catch estimate
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "month_area_recCatch.rds")) %>% 
  # identify months to exclude based on minimal catch or comp estimates 
  #(make region specific)
  mutate(
    min_m = case_when(
      region %in% c("NSoG", "SSoG") ~ 5,
      region %in% c("JdFS", "QCaJS") ~ 6
    ),
    max_m = case_when(
      region == "QCaJS" ~ 8, 
      region %in% c("JdFS", "NSoG", "SSoG")  ~ 9
    ),
    region = fct_relevel(region, "QCaJS", "JdFS", "NSoG", "SSoG")
  ) %>% 
  # drop areas with fewer than 10 datapoints (month-year-area observation = 1)
  group_by(area_n) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(!n < 10) %>%
  droplevels() %>% 
  select(region, region_c, area, area_n, month, month_n, year, catch, eff, 
         eff_z, offset) 


# CLEAN GENETICS  --------------------------------------------------------------

# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                           "recIndProbsLong.rds")) %>% 
  filter(legal == "legal",
         !region %in% c("Queen Charlotte Sound")) %>%
  mutate(
    temp_strata = paste(month, region, sep = "_"),
    sample_id = paste(temp_strata, jDay, year, sep = "_"),
    min_m = case_when(
      region %in% c("N. Strait of Georgia", "S. Strait of Georgia") ~ 5,
      region %in% c("Juan de Fuca Strait",
                    "Queen Charlotte and\nJohnstone Straits") ~ 6
    ),
    max_m = case_when(
      region %in% c("Juan de Fuca Strait", "N. Strait of Georgia",
                    "S. Strait of Georgia") ~ 9,
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
    # ungroup() %>% 
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

# saveRDS(full_dat, here::here("data-raw", "combined_model_inputs.RDS"))


## PREP MODEL INPUTS -----------------------------------------------------------

# add tmb data assuming average effort predictions
# note that these outputs are generated automatically in fit_stockseason
# but aren't retained and are necessary for the current suite of plotting 
# functions
tmb_list <- purrr::pmap(list(catch_dat = full_dat$catch_data, 
                             comp_dat = full_dat$comp_long), 
                        .f = stockseasonr::gen_tmb,
                        random_walk = TRUE,
                        model_type = "integrated")

# join all data into a tibble then generate model inputs
dat <- full_dat %>% 
  mutate(
    comp_wide = purrr::map(tmb_list, "comp_wide"),
    tmb_data = purrr::map(tmb_list, "data"),
    tmb_pars = purrr::map(tmb_list, "pars"),
    tmb_map = purrr::map(tmb_list, "tmb_map"),
    pred_dat_catch = purrr::map(tmb_list, ~ .$pred_dat_catch),
    pred_dat_comp = purrr::map(tmb_list, ~ .$pred_dat_comp)
  ) 

## FIT MODELS ------------------------------------------------------------------

#use fix optim inits except rec_pst (convergence issues)
optim_fix_inits_vec = c(TRUE, FALSE, TRUE, TRUE)

dat2 <- dat %>%
  mutate(sdr = purrr::pmap(list(catch_dat = full_dat$catch_data, 
                                comp_dat = full_dat$comp_long,
                                optim_fix_inits = optim_fix_inits_vec), 
                           .f = stockseasonr::fit_stockseason, 
                           random_walk = TRUE,
                           model_type = "integrated",
                           # nlminb_loops = 2, 
                           silent = FALSE))
dat2$ssdr <- purrr::map(dat2$sdr, summary)

saveRDS(dat2 %>% select(-sdr),
        here::here("generated_data", "model_fits", "combined_model_dir.RDS"))

dat2 <- readRDS(here::here("generated_data", "model_fits",
                           "combined_model_dir.RDS"))


## GENERATE OBS AND PREDICTIONS ------------------------------------------------

# source file for cleaning and plotting functions
source(here::here("R", "functions", "plot_cleaning_functions.R"))
source(here::here("R", "functions", "plot_predictions_splines.R"))

pred_dat <- dat2 %>% 
  mutate(
    #predictions assuming mean effort
    abund_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr), 
                                          .f = gen_abund_pred)),
    cpue_pred_ci = pmap(list(pred_dat_catch, comp_long, ssdr), 
                        gen_abund_pred_area),
    comp_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr), 
                        .f = gen_comp_pred)),
    raw_prop = purrr::map2(comp_long, comp_wide, make_raw_prop_dat),
    raw_abund_area = pmap(list(catch_data, raw_prop), .f = make_raw_abund_dat,
                          spatial_scale = "area"),
    raw_abund_reg = pmap(list(catch_data, raw_prop), .f = make_raw_abund_dat,
                         spatial_scale = "region")
    ) %>% 
  select(dataset:catch_data, abund_pred_ci:raw_abund_reg) %>% 
  mutate(
    comp_pred_ci = map2(grouping_col, comp_pred_ci, .f = stock_reorder),
    raw_prop = map2(grouping_col, raw_prop, .f = stock_reorder)
  )
saveRDS(pred_dat, here::here("generated_data", 
                             "combined_model_predictions.RDS"))
pred_dat <- readRDS(here::here("generated_data",
                               "combined_model_predictions.RDS"))

## PLOT PREDICTIONS ------------------------------------------------------------

#color palette
pal <- readRDS(here::here("generated_data", "disco_color_pal.RDS"))

#pfma map (used to steal legend)
pfma_map <- readRDS(here::here("generated_data", "pfma_map.rds"))

# generic settings
file_path <- here::here("figs", "model_pred", "combined")


## Plot abundance
comm1 <- pred_dat %>% 
  filter(dataset == "gsi_troll_pst") 
comm_abund <- plot_abund(comm1$abund_pred_ci[[1]], 
                         ylab = "Commercial Catch Index") +
  annotate("text", x = -Inf, y = Inf, label = "a)", hjust = -1, vjust = 2)
rec1 <- pred_dat %>% 
  filter(dataset == "gsi_sport_pst") 
rec_abund <- plot_abund(rec1$abund_pred_ci[[1]], 
                        ylab = "Recreational Catch Index") +
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


## Plot area abundance with observations
# combine similar groupings to plot together
comm_area_cpue <- plot_cpue_area(pred_dat$cpue_pred_ci[[1]], 
                        pred_dat$raw_abund_area[[1]])
rec_area_cpue <- plot_cpue_area(pred_dat$cpue_pred_ci[[2]], 
                       pred_dat$raw_abund_area[[2]])


## Composition spline prediction
comp_plots <- map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp, 
                   raw = TRUE)
comp_plots_fix <- map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp, 
                       raw = TRUE, facet_scales = "fixed")

pdf(paste(file_path, "composition_splines.pdf", sep = "/"))
comp_plots
comp_plots_fix
dev.off()


## Composition area prediction
combine_data <- function(in_col = c("pst_agg", "reg1")) {
  pred_dat %>% 
    filter(grouping_col == in_col) %>% 
    select(dataset, grouping_col, comp_pred_ci) %>% 
    unnest(., cols = c(comp_pred_ci))
}  
pst_comp_pred <- combine_data(in_col = "pst_agg")
plot_comp_stacked(pst_comp_pred, grouping_col = "pst_agg")

comp_area_plots <- map2(pred_dat$comp_pred_ci, pred_dat$grouping_col, 
                        plot_comp_stacked)
pdf(paste(file_path, "composition_stacked.pdf", sep = "/"))
comp_area_plots
dev.off()


## Stock-specific catch predictions
ss_abund_plots_free <- map2(pred_dat$comp_pred_ci, pred_dat$raw_abund, 
                            plot_ss_abund, raw = FALSE)
ss_abund_plots_fix <- map2(pred_dat$comp_pred_ci, pred_dat$raw_abund, 
                           plot_ss_abund, raw = FALSE, facet_scales = "fixed")

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
  labs(y = "Predicted Commercial Catch Index")
dev.off()

png(here::here("figs", "ms_figs", "abund_pst_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
ss_abund_plots_free[[2]] +
  labs(y = "Predicted Recreational Catch Index")
dev.off()

png(here::here("figs", "ms_figs", "abund_can_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7)
ss_abund_plots_free[[3]] +
  labs(y = "Predicted Commercial Catch Index")
dev.off()

png(here::here("figs", "ms_figs", "abund_can_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7)
ss_abund_plots_free[[4]] +
  labs(y = "Predicted Recreational Catch Index")
dev.off()

# supplemental figures
png(here::here("figs", "ms_figs", "comm_cpue_fits.png"), res = 400, units = "in",
    height = 5.5, width = 7)
comm_area_cpue
dev.off()

png(here::here("figs", "ms_figs", "rec_cpue_fits.png"), res = 400, units = "in",
    height = 5.5, width = 7)
rec_area_cpue
dev.off()

## Random summary data ---------------------------------------------------------

# Import composition data from fit_dirichlet_splines since extra samples used 
# in recreational analysis
comp_only <- readRDS(here::here("generated_data", "composition_model_data.RDS"))

# How many samples total?
f <- function(x) {
  x %>% 
    select(month, year, region, nn) %>% 
    distinct() %>% 
    pull(nn) %>% 
    sum()
}
f(comp_only$comp_long[[1]])
f(comp_only$comp_long[[2]])

### Summary supplemental tables
tally_f <- function(dat, type, fishery) {
  if (type == "catch") {
    out <- dat %>% 
      group_by(region_c, area, month, month_n) %>% 
      tally() %>% 
      rename(PFMA = area) 
  }
  if (type == "comp") {
    out <- dat %>% 
      group_by(region, month, month_n) %>% 
      tally() 
  }
  out %>% 
    mutate(Fishery = fishery) %>% 
    ungroup()
}

pivot_f <- function(dat1, dat2, type) {
  if (type == "catch") {
    out <- rbind(dat1, dat2) %>% 
      mutate(month = fct_reorder(month, month_n)) %>% 
      arrange(Fishery, region_c, PFMA, month_n) %>% 
      select(-month_n) %>% 
      pivot_wider(names_from = "month", values_from = "n") %>% 
      select(Fishery, Region = region_c, PFMA, `1`:`12`)
  }
  if (type == "comp") {
    out <- rbind(dat1, dat2) %>% 
      mutate(month = fct_reorder(month, month_n)) %>% 
      arrange(Fishery, region, month_n) %>% 
      select(-month_n) %>% 
      pivot_wider(names_from = "month", values_from = "n") %>% 
      select(Fishery, Region = region, `1`:`12`)
  }
  out 
}

# catch data (one observation equals a year-pfma-monthly catch estimate)
rec_n_long <- rec_catch %>% 
  tally_f(., type = "catch", fishery = "rec")
comm_n_long <- comm_catch %>% 
  tally_f(., type = "catch", fishery = "commercial")
catch_n <- pivot_f(rec_n_long, comm_n_long, type = "catch")

# composition data
rec_n2_long <- comp_only$raw_data[[2]] %>% 
  select(id, region, month, month_n) %>% 
  distinct() %>% 
  tally_f(., type = "comp", fishery = "rec")
comm_n2_long <- comp_only$raw_data[[1]] %>% 
  select(id, region, month, month_n) %>% 
  distinct() %>% 
  tally_f(., type = "comp", fishery = "comm")
comp_n <- pivot_f(rec_n2_long, comm_n2_long, type = "comp")

# export
write.csv(catch_n, here::here("figs", "ms_figs", "tableS1_catch_samples.csv"),
          row.names = FALSE)
write.csv(comp_n, here::here("figs", "ms_figs", "tableS2_comp_samples.csv"),
          row.names = FALSE)

# CPUE figs
plot_cpue <- function(dat) {
  dat %>% 
    mutate(cpue = catch / eff,
           log_cpue = log(catch) / log(eff),
           log_cpue2 = log(cpue)) %>% 
    ggplot(., aes(x = month, y = log_cpue2, 
                  fill = fct_reorder(region_c, as.numeric(region)))) +
    geom_point(shape = 21) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    facet_wrap(~fct_reorder(area, as.numeric(region))) +
    ggsidekick::theme_sleek() +
    labs(x = "Month", 
         y = "Catch Per Unit Effort (pieces per boat day)")
}


png(here::here("figs", "ms_figs", "rec_cpue.png"), res = 400, units = "in",
    height = 5.5, width = 7)
plot_cpue(rec_catch)
dev.off()

png(here::here("figs", "ms_figs", "comm_cpue.png"), res = 400, units = "in",
    height = 5.5, width = 7)
plot_cpue(comm_catch)
dev.off()
