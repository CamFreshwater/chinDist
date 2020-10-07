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
make_raw_abund_dat <- function(catch, raw_prop, fishery) {
  if (fishery == "sport") {
    dum <-  catch %>% 
      group_by(region, region_c, month_n, year) %>%
      summarize(
        #sum rec because subareas each have one monthly value
        #mean commercial because daily values for each area
        cum_catch = as.numeric(sum(catch)),
        cum_effort = as.numeric(sum(eff)),
        .groups = "drop") 
  }
  if (fishery == "troll") {
    dum <- catch %>% 
      group_by(region, region_c, area, month_n, year) %>% 
      # calculate daily catch within month first
      summarize(
        mu_catch = mean(catch),
        mu_effort = mean(eff),
        .groups = "drop"
      ) %>% 
      group_by(region, region_c, month_n, year) %>% 
      summarize(
        cum_catch = as.numeric(sum(mu_catch)),
        cum_effort = as.numeric(sum(mu_effort)),
        .groups = "drop"
      ) 
  }
  dum  %>%
    ungroup() %>% 
    mutate(agg_cpue = cum_catch / cum_effort,
           month = as.character(month_n)) %>%
    left_join(., raw_prop, by = c("region", "month_n", "year")) %>%
    mutate(catch_g = samp_g_ppn * cum_catch,
           cpue_g = samp_g_ppn * agg_cpue,
           region = as.factor(region)) %>%
    filter(!is.na(stock))
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
# ------------------------------------------------------------------------------


## Function to generate composition and stock-specific abundance predictions
gen_comp_pred <- function(comp_long, pred_dat_comp, ssdr) {
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
#-------------------------------------------------------------------------------


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
