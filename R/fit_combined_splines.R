## Combined model fits
# Sep 25, 2020
# Fit combined dirichlet/nb model to stock composition and abundance data
# aggregate probability reflects summed probabilities of a given region of 
# origin for a given individual
# Constrain to relatively few months due to limited availability of rec C/E 
# data
# Modified from fit_combined.R to include splines
# Most recent update: Dec. 5 2020

library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)


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
# optim_fix_inits_vec = c(TRUE, FALSE, TRUE, TRUE)
# 
# dat2 <- dat %>%
#   mutate(sdr = purrr::pmap(list(catch_dat = full_dat$catch_data,
#                                 comp_dat = full_dat$comp_long,
#                                 optim_fix_inits = optim_fix_inits_vec),
#                            .f = stockseasonr::fit_stockseason,
#                            random_walk = TRUE,
#                            model_type = "integrated",
#                            # nlminb_loops = 2,
#                            silent = FALSE))
# dat2$ssdr <- purrr::map(dat2$sdr, summary)
# 
# saveRDS(dat2 %>% select(-sdr),
#         here::here("generated_data", "model_fits", "combined_model_dir.RDS"))

dat2 <- readRDS(here::here("generated_data", "model_fits",
                           "combined_model_dir.RDS"))


## GENERATE OBS AND PREDICTIONS ------------------------------------------------

source(here::here("R", "functions", "plot_cleaning_functions.R"))

pred_dat <- dat2 %>%
  mutate(
    #predictions assuming mean effort
    abund_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr),
                                          .f = gen_abund_pred)),
    cpue_pred_ci = pmap(list(pred_dat_catch, comp_long, ssdr),
                        gen_abund_pred_area),
    abund_year_ci = pmap(list(ssdr, catch_data, pred_dat_catch, cpue_pred_ci),
                         gen_rand_int_pred),
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

source(here::here("R", "functions", "plot_predictions_splines.R"))

#color palette
pal <- readRDS(here::here("generated_data", "disco_color_pal.RDS"))

#pfma map (used to steal legend)
pfma_map <- readRDS(here::here("generated_data", "pfma_map.rds"))

# generic settings
file_path <- here::here("figs", "model_pred", "combined")


## Plot abundance - fixed effects only
comm1 <- pred_dat %>% 
  filter(dataset == "gsi_troll_pst") 
comm_abund <- plot_abund(comm1$abund_pred_ci[[1]], 
                         ylab = "Commercial Predicted\nStandardized CPUE") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -1, vjust = 2)
rec1 <- pred_dat %>% 
  filter(dataset == "gsi_sport_pst") 
rec_abund <- plot_abund(rec1$abund_pred_ci[[1]], 
                        ylab = "Recreational Predicted\nStandardized CPUE") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -1, vjust = 2)

combo_abund <- cowplot::plot_grid(comm_abund, rec_abund, nrow = 2, 
                                  align = "v") %>% 
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


## Plot abundance - include random effects
# combine data 
comm_yr_catch <- comm1$abund_year_ci[[1]] %>% 
  group_by(month_n, region_c, year) %>% 
  summarize(pred_catch_yr = sum(pred_catch_yr), .groups = "drop")  %>% 
  plot_catch_yr(., facet_scale = "region_c") 
rec_yr_catch <- rec1$abund_year_ci[[1]] %>% 
  group_by(month_n, region_c, year) %>% 
  summarize(pred_catch_yr = sum(pred_catch_yr), .groups = "drop") %>% 
  plot_catch_yr(., facet_scale = "region_c")


## Plot area abundance with observations
# combine similar groupings to plot together
comm_area_cpue <- plot_cpue_area(pred_dat$cpue_pred_ci[[1]], 
                        pred_dat$raw_abund_area[[1]])
rec_area_cpue <- plot_cpue_area(pred_dat$cpue_pred_ci[[2]], 
                       pred_dat$raw_abund_area[[2]])


## Composition spline prediction (now occurs in fit_dirichlet_splines to maximize
# monthly coverage)
# comp_plots <- map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp,
#                    raw = TRUE)
# comp_plots_fix <- map2(pred_dat$comp_pred_ci, pred_dat$raw_prop, plot_comp, 
#                        raw = TRUE, facet_scales = "fixed")
# 
# pdf(paste(file_path, "composition_splines.pdf", sep = "/"))
# comp_plots
# comp_plots_fix
# dev.off()


## Composition area prediction (now occurs in fit_dirichlet_splines to maximize
# monthly coverage)
# combine_data <- function(in_col = c("pst_agg", "reg1")) {
#   pred_dat %>%
#     filter(grouping_col == in_col) %>%
#     select(dataset, grouping_col, comp_pred_ci) %>%
#     unnest(., cols = c(comp_pred_ci))
# }
# pst_comp_pred <- combine_data(in_col = "pst_agg")
# pst_area <- plot_comp_stacked(pst_comp_pred, grouping_col = "pst_agg", 
#                               palette_name = "rainbow")
# can_comp_pred <- combine_data(in_col = "reg1")
# can_area <- plot_comp_stacked(can_comp_pred, grouping_col = "reg1")
# 
# pdf(paste(file_path, "composition_stacked.pdf", sep = "/"))
# pst_area
# can_area
# dev.off()


## Stock-specific catch predictions
ss_abund_plots_free <- map2(pred_dat$comp_pred_ci, pred_dat$raw_abund_reg, 
                            plot_ss_abund, raw = FALSE)
ss_abund_plots_fix <- map2(pred_dat$comp_pred_ci, pred_dat$raw_abund_reg, 
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

## REPLACED WITH ESTIMATES FROM COMPOSITION ONLY MODEL
# png(here::here("figs", "ms_figs", "comp_pst_comm.png"), res = 400, units = "in",
#     height = 5.5, width = 7.5)
# comp_plots_fix[[1]]
# dev.off()
# 
# png(here::here("figs", "ms_figs", "comp_pst_rec.png"), res = 400, units = "in",
#     height = 5.5, width = 7.5)
# comp_plots_fix[[2]]
# dev.off()
# 
# png(here::here("figs", "ms_figs", "comp_can_comm.png"), res = 400, units = "in",
#     height = 5.5, width = 7)
# comp_plots_fix[[3]]
# dev.off()
# 
# png(here::here("figs", "ms_figs", "comp_can_rec.png"), res = 400, units = "in",
#     height = 5.5, width = 7)
# comp_plots_fix[[4]]
# dev.off()

png(here::here("figs", "ms_figs", "abund_pst_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
ss_abund_plots_free[[1]] +
  labs(y = "Commercial Predicted Standardized CPUE")
dev.off()

png(here::here("figs", "ms_figs", "abund_pst_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7.5)
ss_abund_plots_free[[2]] +
  labs(y = "Recreational Predicted Standardized CPUE")
dev.off()

png(here::here("figs", "ms_figs", "abund_can_comm.png"), res = 400, units = "in",
    height = 5.5, width = 7)
ss_abund_plots_free[[3]] +
  labs(y = "Commercial Predicted Standardized CPUE")
dev.off()

png(here::here("figs", "ms_figs", "abund_can_rec.png"), res = 400, units = "in",
    height = 5.5, width = 7)
ss_abund_plots_free[[4]] +
  labs(y = "Recreational Predicted Standardized CPUE")
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

png(here::here("figs", "ms_figs", "comm_yr_catch.png"), res = 400, units = "in",
    height = 5.5, width = 7)
comm_yr_catch
dev.off()

png(here::here("figs", "ms_figs", "rec_yr_catch.png"), res = 400, units = "in",
    height = 5.5, width = 7)
rec_yr_catch
dev.off()



## Supplementary summary data --------------------------------------------------

# Import composition data from fit_dirichlet_splines since extra samples used 
# in recreational analysis
comp_only <- readRDS(here::here("generated_data", 
                                "composition_model_predictions.RDS"))

# How many samples total?
f <- function(x) {
  x %>% 
    select(month, year, region, nn) %>% 
    distinct() %>% 
    pull(nn) %>% 
    sum()
}

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

# catch data (one observation equals a year-pfma-monthly catch estimate)
rec_n_long <- rec_catch %>% 
  tally_f(., type = "catch", fishery = "rec")
comm_n_long <- comm_catch %>% 
  tally_f(., type = "catch", fishery = "commercial")

catch_n_list <- rbind(rec_n_long, comm_n_long) %>% 
  split(., .$region_c) 
plot_f <- function(dat) {
  alpha_labs <- unique(dat$PFMA)
  alpha_vals <- seq(0.4, 1, length = length(alpha_labs))
  ggplot(dat) +
    geom_bar(aes(x = month_n, y = n, fill = region_c, alpha = PFMA), 
             position = position_dodge(0.9),
             stat = "identity") +
    scale_fill_manual(name = "Region", values = pal, guide = FALSE) +
    scale_alpha_manual(labels = alpha_labs, values = alpha_vals) +
    scale_x_continuous(breaks = seq(2, 12, by = 2), limits = c(1, 12)) +
    scale_y_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 11)) +
    ggsidekick::theme_sleek() +
    theme(plot.margin=unit(c(1.5, 1.5, 1.5, 1.5), "points"),
          axis.title = element_blank())
}

catch_n_plots <- map(catch_n_list, plot_f)
catch_n_plot <- cowplot::plot_grid(
  catch_n_plots[["NWVI"]], catch_n_plots[["SWVI"]],
  catch_n_plots[["Queen Charlotte and\nJohnstone Straits"]], 
  catch_n_plots[["N. Strait of Georgia"]], 
  catch_n_plots[["S. Strait of Georgia"]],
  catch_n_plots[["Juan de Fuca Strait"]],
  nrow = 3, align = "v"
) %>% 
  arrangeGrob(., 
              bottom = textGrob("Month", 
                                gp = gpar(col = "grey30", fontsize=11)),
              left = textGrob("Years with Samples", 
                              gp = gpar(col = "grey30", fontsize=11), 
                              rot = 90)) %>% 
  grid.arrange()
# add legend from PFMA map
catch_n_plot2 <- cowplot::plot_grid(
  cowplot::get_legend(pfma_map),
  catch_n_plot,
  ncol=1, rel_heights=c(.1, .9)
)

png(here::here("figs", "ms_figs", "catch_sample_sizes.png"), res = 400, 
    units = "in",
    height = 5.5, width = 7)
catch_n_plot2
dev.off()


# composition data
rec_n2_long <- comp_only$raw_data[[2]] %>% 
  select(id, region, month, month_n) %>% 
  distinct() %>% 
  tally_f(., type = "comp", fishery = "rec")
comm_n2_long <- comp_only$raw_data[[1]] %>% 
  select(id, region, month, month_n) %>% 
  distinct() %>% 
  tally_f(., type = "comp", fishery = "comm")
comp_n_dat <- rbind(rec_n2_long, comm_n2_long) %>% 
  mutate(region = fct_relevel(region, "NWVI", "SWVI", 
                              "Queen Charlotte and\nJohnstone Straits",
                              "Juan de Fuca Strait", "N. Strait of Georgia",
                              "S. Strait of Georgia"))
comp_n_plot <- ggplot(comp_n_dat) +
  geom_bar(aes(x = month, y = n, fill = region), 
           stat = "identity") +
  scale_fill_manual(name = "Region", values = pal) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top") +
  labs(y = "Number of Sampled Individuals", x = "Month") +
  facet_wrap(~region)

png(here::here("figs", "ms_figs", "comp_sample_sizes.png"), res = 400, 
    units = "in",
    height = 5.5, width = 7)
comp_n_plot
dev.off()


# REPLACED BY ANALYSIS IN MANAGEMENT_RESTRICTIONS_SUPP>.Rmd
# # composition samples by year, divided among kept and retained
# tally_gsi <- function(dat) {
#   gear <- dat$gear
#   if (any(gear == "troll")) {
#     out <- dat %>% 
#       mutate(release = "Kept")
#   } else {
#     out <- dat
#   }
#   out %>%
#     group_by(gear, region, month_n, year, release) %>% 
#     summarize(count = sum(adj_prob), .groups = "drop")
# } 
# 
# tally_dat <- map(comp_only$raw_data, tally_gsi) %>% 
#   bind_rows()  %>%
#   filter(gear == "sport",
#          #constrain to months when management measures to protect spring
#          #runs are in place
#          month_n > 3, 
#          month_n < 8) %>% 
#   group_by(region, month_n, release, year) %>% 
#   summarize(sampled_ind = sum(count), .groups = "drop") %>% 
#   mutate(release = fct_recode(release, Released = "Rel"),
#          region = fct_relevel(region, 
#                               "Queen Charlotte and\nJohnstone Straits", 
#                               after = 0))
# 
# kept_released_tally <- expand.grid(region = unique(tally_dat$region),
#                                    year = unique(tally_dat$year),
#                                    release = unique(tally_dat$release),
#                                    month_n = unique(tally_dat$month_n)) %>% 
#   arrange(region, year, release, month_n) %>% 
#   left_join(., tally_dat, by = c("region", "year", "release", "month_n")) %>% 
#   mutate(month = as.factor(month_n))
# 
# kept_released <- ggplot(kept_released_tally) +
#   geom_bar(aes(x = year, y = sampled_ind, fill = release), 
#                position = position_dodge(0.9),
#            stat = "identity"
#            ) +
#   labs(x = "Year", y = "Number of Individuals") +
#   scale_fill_manual(values = c("#999999", "#E69F00"),
#                     name = "") +
#   facet_wrap(~region) +
#   ggsidekick::theme_sleek()
# 
# png(here::here("figs", "ms_figs", "kept_retained_samples.png"), res = 400, 
#     units = "in",
#     height = 5.5, width = 7)
# kept_released
# dev.off()
