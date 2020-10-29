## Miscellaneous cleaning functions to generate predictions and prep data for 
# plotting.
# Only tested with fit_combined_splines and fit_dirichlet_splines
# Sep. 26, 2020

## Function to make raw proportion data
make_raw_prop_dat <- function(comp_long, comp_wide) {
  n_stks <- length(unique(comp_long$agg))
  comp_wide %>%
    pivot_longer((ncol(.) - n_stks + 1):ncol(.), 
                 names_to = "agg", values_to = "agg_prob") %>% 
    group_by(region, region_c, month_n, year, agg) %>%
    summarize(samp_g = sum(agg_prob),
              .groups = "drop") %>%
    group_by(region, region_c, month_n, year) %>%
    mutate(samp_total = sum(samp_g)) %>% 
    ungroup() %>% 
    mutate(samp_g_ppn = samp_g / samp_total,
           stock = fct_reorder(agg, desc(samp_g_ppn))) 
}
# ------------------------------------------------------------------------------

## Function to calculate raw stock-specific abundance observations
# data_type necessary for monthly vs. daily predictions
make_raw_abund_dat <- function(catch, raw_prop, 
                               spatial_scale = c("area", "region")) {
  if (spatial_scale == "area") {
    out <- catch %>% 
      mutate(area_cpue = catch / eff,
             log_area_cpue = log(catch) / offset)  
  }
  if (spatial_scale == "region") {
    out <- catch %>% 
      group_by(region, region_c, month_n, year) %>%
      summarize(reg_catch = sum(catch),
                reg_eff = sum(eff),
                log_cpue = log(reg_catch) / log(reg_eff),
                month = as.character(month_n),
                .groups = "drop") %>%
      ungroup() %>%
      distinct() %>% 
      left_join(., raw_prop, by = c("region", "region_c", "month_n", "year")) %>%
      mutate(log_cpue_g = samp_g_ppn * log_cpue,
             region = as.factor(region)) %>%
      filter(!is.na(stock)) 
  }
  return(out)
}
# ------------------------------------------------------------------------------

## Function to generate aggregate abundance predictions (use comp predictions
# because aggregate)
gen_abund_pred <- function(comp_long, pred_dat_comp, ssdr) {
  abund_pred <- ssdr[rownames(ssdr) %in% "log_agg_pred_abund", ] 
  pred_dat <- pred_dat_comp %>% 
    left_join(., comp_long %>% select(region, region_c) %>% distinct(), 
              by = "region") %>% 
    mutate(region_c = fct_reorder(region_c, as.numeric(region)))
  
  data.frame(est_link = abund_pred[ , "Estimate"],
             se_link =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_dat, .) %>% 
    mutate(
      pred_est = exp(est_link),
      pred_low = exp(est_link + (qnorm(0.025) * se_link)),
      pred_up = exp(est_link + (qnorm(0.975) * se_link))
    )
}

# like above but for areas (hence pred_dat_catch) and includes cpue instead of 
# stock-specific est
gen_abund_pred_area <- function(pred_dat_catch, comp_long, ssdr) {
  abund_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] 
  pred_dat <- pred_dat_catch %>% 
    left_join(., comp_long %>% select(region, region_c) %>% distinct(), 
              by = "region") %>% 
    mutate(region_c = fct_reorder(region_c, as.numeric(region)),
           offset = unique(pred_dat_catch$offset))
  
  data.frame(est_link = abund_pred[ , "Estimate"],
             se_link =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_dat, .) %>% 
    mutate(
      pred_est = exp(est_link),
      pred_low = exp(est_link + (qnorm(0.025) * se_link)),
      pred_up = exp(est_link + (qnorm(0.975) * se_link)),
      pred_est_logcpue = est_link / offset,
      pred_low_logcpue = (est_link + (qnorm(0.025) * se_link)) / offset,
      pred_up_logcpue = (est_link + (qnorm(0.975) * se_link)) / offset
    )
}
# ------------------------------------------------------------------------------

## Function to incorporate random intercepts in area-specific abundance 
# estimates
gen_rand_int_pred <- function(ssdr, catch_dat, pred_dat_catch, cpue_pred) {
  # save random intercept parameters and bind to vector of years based on input
  # data
  rand_int <- ssdr[rownames(ssdr) %in% "z1_k", ] 
  year_ints <- data.frame(year = levels(as.factor(as.character(catch_dat$year))),
                          z1_k = as.numeric(rand_int[, "Estimate"])) 
  # generate predictions
  year_preds <- expand.grid(year = year_ints$year, 
                            month = unique(pred_dat_catch$month)) %>% 
    left_join(pred_dat_catch %>% select(area, month), ., by = "month") %>% 
    left_join(., year_ints, by = "year")
  
  cpue_pred %>% 
    select(month_n:est_link) %>% 
    left_join(., year_preds, by = c("area", "month")) %>% 
    mutate(
      year = as.factor(year),
      est_link_yr = est_link + z1_k,
      pred_est_logcpue = est_link_yr / offset,
      pred_catch_yr = exp(est_link_yr)
    ) 
}


#-------------------------------------------------------------------------------

pred_dat_comp <- dat2$pred_dat_comp[[1]]
comp_long <- dat2$comp_long[[1]]
ssdr <- dat2$ssdr[[1]]

## Function to generate composition and stock-specific abundance predictions
gen_comp_pred <- function(comp_long, pred_dat_comp, ssdr, comp_only = FALSE) {
  comp_pred <- ssdr[rownames(ssdr) %in% "inv_logit_pred_pi_prop", ]
  comp_abund_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund_mg", ]
  beta_comp <- ssdr[rownames(ssdr) %in% "b2_jg", ]
  
  stk_names <- unique(comp_long$agg)
  n_preds <- nrow(pred_dat_comp)
  
  pred_dat <- pred_dat_comp %>% 
    left_join(., comp_long %>% select(region, region_c) %>% distinct(), 
              by = "region") %>%
    mutate(region_c = fct_reorder(region_c, as.numeric(region)))
  dum <- purrr::map_dfr(seq_along(stk_names), ~ pred_dat)
  
  if (comp_only) {
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
  } else {
    data.frame(
      stock = as.character(rep(stk_names, each = n_preds)),
      link_prob_est = comp_pred[ , "Estimate"],
      link_prob_se =  comp_pred[ , "Std. Error"],
      link_abund_est = comp_abund_pred[ , "Estimate"],
      link_abund_se =  comp_abund_pred[ , "Std. Error"]
    ) %>% 
      cbind(dum, .) %>%
      mutate(
        pred_prob_est = car::logit(link_prob_est),
        pred_prob_low = pmax(0,
                             car::logit(link_prob_est + (qnorm(0.025) *
                                                           link_prob_se))),
        pred_prob_up = car::logit(link_prob_est + (qnorm(0.975) * link_prob_se)),
        comp_abund_est = exp(link_abund_est) / 1000,
        comp_abund_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)) / 1000,
        comp_abund_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se)) / 1000
      ) 
  }
  
}
#-------------------------------------------------------------------------------

## Function to reorder stocks for plotting
stock_reorder <- function(key, comp_data) {
  if (key == "pst_agg") {
    out <- comp_data %>% 
      mutate(stock =  fct_relevel(stock, "NBC_SEAK", "WCVI", "FR-early", 
                                  "FR-late", "SOG", "PSD", "WACST", "CR-sp&su", 
                                  "CR-bright", "CR-tule", "CA_ORCST"),
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
