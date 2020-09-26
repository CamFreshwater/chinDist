## Dirichlet model fit
# Sep 25, 2020
# Fit combined dirichlet/nb model to stock composition data
# Duplicates composition component of combined model. Ran independently to 
# increase temporal range
# Same as fit_dirichlet but adds splines 

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

tt <- rec %>%
  select(region, id, month_n, year) %>%
  distinct()
table(tt$region, tt$month_n)

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
      mutate(region = fct_relevel(region, "QCaJS", "NSoG", "SSoG", "JdFS"))
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

# comp_in <- comp$data[[2]]
# data_type <- comp$dataset[[2]]

prep_dir_inputs <- function(comp_in, data_type) {
  gsi_trim <- comp_in %>% 
    group_by(region, month, year) %>%
    mutate(nn = sum(agg_prob)) %>%
    #remove strata with less than 10 individuals total 
    filter(!nn < 10) %>%
    ungroup() %>%
    droplevels() %>%
    select(sample_id, region, region_c, year, month, month_n, agg, agg_prob) %>%
    distinct()
  
  gsi_wide <- gsi_trim %>% 
    pivot_wider(., names_from = agg, values_from = agg_prob) %>%
    mutate_if(is.numeric, ~replace_na(., 0.000001))
  
  y_obs <- gsi_wide %>% 
    select(-c(sample_id:month_n)) %>% 
    as.matrix() 

  yr_vec <- as.numeric(gsi_wide$year) - 1
  
  #generate model matrix based on GAM
  months <- unique(gsi_wide$month_n)
  n_months <- length(months)
  n_knots <- ifelse(max(months) > 8, 4, 3)
  spline_type <- ifelse(max(months) == 12, "cc", "tp")
  # response variable doesn't matter, since not fit
  m1 <- gam(rep(0, length.out = nrow(gsi_wide)) ~ 
              region + s(month_n, bs = spline_type, k = n_knots, by = region),
            data = gsi_wide)
  fix_mm <- predict(m1, type = "lpmatrix")
  
  # data frame for predictions
  # account for strong differences in sampling months for sport fishery
  if (comp_in$gear[1] == "sport") {
    pred_dat <- split(gsi_wide, gsi_wide$region) %>% 
      map(., function (x) {
        expand.grid(
          month_n = seq(min(x$month_n), 
                        max(x$month_n),
                        by = 0.1),
          region = unique(x$region)
        )
      }) %>% 
      bind_rows()
  } else {
    pred_dat <- expand.grid(
      month_n = seq(min(gsi_wide$month_n), 
                    max(gsi_wide$month_n),
                    by = 0.1),
      region = unique(gsi_wide$region)
    )
  }
  pred_mm <- predict(m1, pred_dat, type = "lpmatrix")
  
  data <- list(y_obs = y_obs, #obs
               rfac = yr_vec, #random intercepts
               fx_cov = fix_mm, #fixed cov model matrix
               n_rfac = length(unique(yr_vec)), #number of random intercepts
               pred_cov = pred_mm
  ) 
  parameters <- list(z_ints = matrix(0, nrow = ncol(fix_mm), 
                                     ncol = ncol(y_obs)),
                     z_rfac = rep(0, times = length(unique(yr_vec))),
                     log_sigma_rfac = 0
                     )
  
  list("fix_mm" = fix_mm, "pred_dat" = pred_dat, 
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "long_data" = gsi_trim, "wide_data" = gsi_wide)
}

# add model parameters
comp2 <- comp %>%
  mutate(model_inputs = map2(comp_long, dataset, .f = prep_dir_inputs),
         # reformat comp2 to match dat2 from fit_combined_splines to share 
         # plotting functions
         comp_wide = map(model_inputs, "wide_data"),
         pred_dat_comp = map(model_inputs, "pred_dat")) 


# FIT --------------------------------------------------------------------------

compile(here::here("src", "dirichlet_randInt.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_randInt")))

fit_model <- function(x) {
  ## Make a function object
  # x <- comp2$model_inputs[[2]]
  obj <- MakeADFun(data = x$data, 
                   parameters = x$parameters, 
                   random = c("z_rfac"),
                   DLL = "dirichlet_randInt"
  )
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  f_name <- paste(x$data_type, "dirichlet_spline_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generated_data", "model_fits", f_name))
}

# fit
map(comp2$model_inputs, .f = fit_model)


## PLOT PREDICTIONS ------------------------------------------------------------

comp2$ssdr <- map(comp2$model_inputs, function (x) {
  f_name <- paste(x$data_type, "dirichlet_spline_ssdr.RDS", sep = "_")
  readRDS(here::here("generated_data", "model_fits", f_name))
})


# source file for cleaning and plotting functions
source(here::here("R", "functions", "plot_cleaning_functions.R"))
source(here::here("R", "functions", "plot_predictions_splines.R"))


# generate predictions for plotting (stripped down version of 
# gen_comp_pred)
gen_comp_only_pred <- function(comp_long, pred_dat_comp, ssdr) {
  comp_pred <- ssdr[rownames(ssdr) %in% "inv_logit_pred_pi_prop", ]
  beta_comp <- ssdr[rownames(ssdr) %in% "z_ints", ]
  
  stk_names <- unique(comp_long$agg)
  n_preds <- nrow(pred_dat_comp)
  
  pred_dat <- pred_dat_comp %>% 
    left_join(., comp_long %>% select(region, region_c) %>% distinct(), 
              by = "region") %>%
    mutate(region_c = fct_reorder(region_c, as.numeric(region)))
  dum <- purrr::map_dfr(seq_along(stk_names), ~ pred_dat)
  
  data.frame(
    stock = as.character(rep(stk_names, each = n_preds)),
    link_prob_est = comp_pred[ , "Estimate"],
    link_prob_se =  comp_pred[ , "Std. Error"]
  ) %>% 
    cbind(dum, .) %>%
    mutate(
      pred_prob_est = car::logit(link_prob_est),
      pred_prob_low = pmax(0,
                           car::logit(link_prob_est + (qnorm(0.025) *
                                                         link_prob_se))),
      pred_prob_up = car::logit(link_prob_est + (qnorm(0.975) * link_prob_se))
    )  
}

pred_dat_comp <- comp2 %>% 
  mutate(
    comp_pred_ci = suppressWarnings(pmap(list(comp_long, pred_dat_comp, ssdr), 
                                         .f = gen_comp_only_pred)),
    raw_prop = purrr::map2(comp_long, comp_wide, make_raw_prop_dat)
    ) %>%
  select(dataset, grouping_col, comp_pred_ci, raw_prop) %>% 
  mutate(
    comp_pred_ci = map2(grouping_col, comp_pred_ci, .f = stock_reorder),
    raw_prop = map2(grouping_col, raw_prop, .f = stock_reorder)
  ) 


#color palette
# pal <- readRDS(here::here("generated_data", "color_pal.RDS"))
pal <- readRDS(here::here("generated_data", "disco_color_pal.RDS"))


comp_plots <- map2(pred_dat_comp$comp_pred_ci, pred_dat_comp$raw_prop, plot_comp, 
                   raw = TRUE)
comp_plots_fix <- map2(pred_dat_comp$comp_pred_ci, pred_dat_comp$raw_prop, 
                       plot_comp, raw = TRUE, facet_scales = "fixed")
comp_plots_fix2 <- map2(pred_dat_comp$comp_pred_ci, pred_dat_comp$raw_prop, 
                        plot_comp, raw = FALSE, facet_scales = "fixed")

pdf(here::here("figs", "model_pred", "dirichlet_only", "composition_splines.pdf"))
comp_plots
dev.off()

pdf(here::here("figs", "model_pred", "dirichlet_only", 
               "composition_splines_fixed.pdf"))
comp_plots_fix
dev.off()


# MS figs
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
