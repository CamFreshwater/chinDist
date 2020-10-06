## Principal Tensor Analysis 
# Oct. 6 2020
# Explore using PTA to identify clusters of stocks with similar characteristics
# using predicted abundance from splines

library(tidyverse)
library(ade4)
library(PTAk)

pred_dat <- readRDS(here::here("generated_data",
                               "combined_model_predictions.RDS"))

# convert predictions to appropriate format
# 1. scale estimates of abundance WITHIN each region to account for differences
# in catchability, then split into a list by region

pred_list <- pred_dat %>% 
  filter(dataset %in% c("gsi_troll_pst", "gsi_sport_pst")) %>% 
  select(dataset, comp_pred_ci) %>% 
  mutate(scaled_pred = map(comp_pred_ci, function (x) {
    x %>% 
      select(month_n, region, stock, comp_abund_est) %>% 
      mutate(scale_abund = scale(comp_abund_est, center = TRUE, 
                                 scale = TRUE)[ , 1]) 
  })) %>% 
  select(-comp_pred_ci) %>% 
  unnest(scaled_pred) %>% 
  select(-dataset, -comp_abund_est) %>%
  # either filter to remove months or filter fewer months and drop some regions
  filter(!month_n < 6,
         !month_n > 8) %>% 
  #scale within a stock to generate anomalies
  group_by(stock) %>% 
  mutate(scale_abund = scale(scale_abund, center = TRUE, scale = TRUE)[ , 1]) %>%
  ungroup() %>% 
  glimpse()
  pivot_wider(names_from = stock, values_from = scale_abund) %>% 
  split(., .$region) 
