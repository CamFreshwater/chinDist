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
    region = fct_relevel(region, "QCaJS", "NSoG", "SSoG", "JdFS")
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
      mutate(region = fct_relevel(region, "QCaJS", "NSoG", "SSoG", "JdFS"))
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

# catch_dat <- full_dat$catch_data[[4]]
# comp_dat <- full_dat$comp_long[[4]]

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
              offset,
            data = catch_dat,
            family = nb)
  fix_mm_catch <- predict(m1, type = "lpmatrix")

  # position of offset vector
  offset_pos <- grep("^offset$", colnames(fix_mm_catch))
  
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
    offset = mean(catch_dat$offset)
  ) %>%
    left_join(., 
              catch_dat %>% select(region, area) %>% distinct(), 
              by = "area") %>% 
    mutate(strata = paste(month_n, region, sep = "_")) %>% 
    filter(strata %in% comp_strata) %>% 
    arrange(region) %>%
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
    mutate_if(is.numeric, ~replace_na(., 0.00001))
  
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
            data = gsi_wide)
  fix_mm_comp <- predict(m2, type = "lpmatrix")
  
  pred_mm_comp <- predict(m2, pred_dat, type = "lpmatrix")

  ## region-stock combinations with 0 observations to map
  comp_map <- comp_dat %>%
    group_by(region, agg) %>%
    complete(., region, nesting(agg)) %>% 
    summarize(total_obs = sum(agg_prob), .groups = "drop") %>% 
    filter(is.na(total_obs)) %>% 
    select(agg, region)
  
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
    b2_jg = matrix(runif(ncol(fix_mm_comp) * ncol(obs_comp), 0.1, 0.9), 
                   nrow = ncol(fix_mm_comp), 
                   ncol = ncol(obs_comp)),
    z2_k = rep(0, times = length(unique(yr_vec_comp))),
    log_sigma_zk2 = log(0.5)
  )
  
  # mapped values (not estimated)
  # offset for effort
  pars$b1_j[offset_pos] <- 1
  b1_j_map <- seq_along(pars$b1_j)
  b1_j_map[offset_pos] <- NA
  tmb_map <- list(b1_j = as.factor(b1_j_map))
  # offset for composition parameterss 
  if (!is.na(comp_map$agg[1])) {
    temp_betas <- pars$b2_jg
    for (i in 1:nrow(comp_map)) {
      offset_stock_pos <- grep(paste(comp_map$agg[1], collapse="|"), 
                               colnames(obs_comp))
      offset_region_pos <- grep(paste(comp_map$region[1], collapse="|"), 
                                colnames(fix_mm_comp))
      for (j in seq_len(ncol(fix_mm_comp))) {
        for (k in seq_len(ncol(obs_comp))) {
          if (j %in% offset_region_pos & k %in% offset_stock_pos) {
            pars$b2_jg[j, k] <- 0.00001
            temp_betas[j, k] <- NA
          }
        }
      }
    }
    comp_map_pos <- which(is.na(as.vector(temp_betas)))
    b2_jg_map <- seq_along(pars$b2_jg)
    b2_jg_map[comp_map_pos] <- NA
    tmb_map <- c(tmb_map, list(b2_jg = as.factor(b2_jg_map)))
  }

  list("data" = data, "pars" = pars, "comp_wide" = gsi_wide,
       "pred_dat_catch" = pred_dat_catch, "pred_dat_comp" = pred_dat,
       "tmb_map" = tmb_map)
}


#add tmb data assuming average effort predictions (for month specific look at
# non-spline code)
tmb_list <- map2(full_dat$catch_data, full_dat$comp_long, .f = gen_tmb)

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

compile(here::here("src", "nb_dirichlet_1re.cpp"))
dyn.load(dynlib(here::here("src", "nb_dirichlet_1re")))

fit_mod <- function(tmb_data, tmb_pars, tmb_map, nlminb_loops = 2) {
  obj <- MakeADFun(tmb_data, tmb_pars, map = tmb_map,
                   #dat$tmb_data[[4]], dat$tmb_pars[[4]], map = dat$tmb_map[[4]],
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
  mutate(sdr = purrr::pmap(list(tmb_data, tmb_pars, tmb_map), 
                           .f = fit_mod, nlminb_loops = 2),
         ssdr = purrr::map(sdr, summary)
         )

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
    comp_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr), 
                        .f = gen_comp_pred)),
    raw_prop = purrr::map2(comp_long, comp_wide, make_raw_prop_dat),
    raw_abund = pmap(list(catch_data, raw_prop), .f = make_raw_abund_dat)
    ) %>% 
  select(dataset:catch_data, abund_pred_ci:raw_abund) %>% 
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
    mutate(cpue = catch / eff) %>% 
    ggplot(., aes(x = month, y = cpue, 
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
