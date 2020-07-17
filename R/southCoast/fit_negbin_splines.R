## Neg bin model fit w/ splines
# July 20, 2020
# Fit neg binomial model to catch data from WCVI troll fishery and southern BC
# rec fisheries
# Modified from fit_negbin.R to include splines

library(tidyverse)
library(TMB)
library(mgcv)

# clean function 
clean_catch <- function(dat) {
  dat %>% 
    filter(!is.na(eff),
           !eff == "0") %>% 
    mutate(region = factor(region),
           area_n = as.numeric(area),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff)),
           eff_z2 = eff_z^2) %>% 
    arrange(region, month) 
}

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  rename(eff = boatDays, region = catchReg) %>% 
  clean_catch(.) %>% 
  # drop inside areas where seasonal catches not available
  filter(!area_n < 100) %>% 
  droplevels()

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
            eff = sum(subarea_eff),
            cpue = catch / eff) %>% 
  ungroup() %>% 
  mutate(temp_strata = paste(month_n, region, sep = "_"),
         region_c = as.character(region),
         region = abbreviate(region, minlength = 4))  %>% 
  clean_catch(.) %>% 
  # drop months with minimal catch estimates
  filter(!month_n < 5,
         !month_n > 9) %>% 
  # drop areas with fewer than 10 datapoints
  group_by(area_n) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(!n < 10) %>% 
  droplevels()

table(rec_catch$month_n, rec_catch$area, rec_catch$region)

# PREP DATA --------------------------------------------------------------------

prep_catch <- function (catch_dat, data_type = NULL, kk = 4) {
  yr_vec <- as.numeric(as.factor(as.character(catch_dat$year))) - 1
  
  # model matrix for fixed effects
  # fix_mm1 <- model.matrix(~ area + month + eff_z + eff_z2, catch_dat)
  
  #generate model matrix based on GAM
  months <- unique(catch_dat$month_n)
  n_months <- length(months)
  spline_type <- ifelse(n_months == 12, "cc", "tp")
  m1 <- gam(catch ~ s(month_n, bs = spline_type, k = kk, by = area) +
                   eff_z + eff_z2, 
                 knots = list(month_n = c(min(months), max(months))),
                 data = catch_dat,
                 family = nb)
  fix_mm <- predict(m1, type = "lpmatrix")
  
  # make predictive model matrix including null values for effort
  pred_dat <- expand.grid(
    month_n = seq(min(catch_dat$month_n), 
                  max(catch_dat$month_n),
                  length.out = 50),
    area = unique(catch_dat$area),
    eff_z = 0,
    eff_z2 = 0
  ) %>% 
    left_join(., catch_dat %>% select(region, area), by = "area") %>% 
    mutate(reg_month = paste(region, round(month_n, 3), sep = "_")) %>% 
    distinct()
  pred_mm_catch <- predict(m1, pred_dat, type = "lpmatrix")
  
  # construct factor key for regional aggregates associated with areas
  # pred_dat_catch <- make_pred_dat(catch_dat) 
  grouping_vec <- as.numeric(as.factor(as.character(pred_dat$reg_month))) - 1
  grouping_key <- unique(grouping_vec)
  
  data <- list(y1_i = catch_dat$catch,
               X1_ij = fix_mm,
               factor1k_i = yr_vec,
               nk1 = length(unique(yr_vec)),
               X1_pred_ij = pred_mm_catch,
               pred_factor2k_h = grouping_vec,
               pred_factor2k_levels = grouping_key
  )
  
  parameters <- list(
    #abundance parameters
    b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
    log_phi = log(1.5),
    z1_k = rep(0, data$nk1),
    log_sigma_zk1 = log(0.25)
  )
  
  if (is.null(data_type)) {
    data_type <- unique(catch_dat$legal)
  }
  
  list("fix_mm" = fix_mm, "mm_pred" = pred_mm_catch, "pred_data" = pred_dat, 
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "input_data" = catch_dat, "m1" = m1)
}

comm_list <- prep_catch(comm_catch, data_type = "comm")
rec_list <- prep_catch(rec_catch, data_type = "rec")
fishery_list <- list(comm_list, rec_list)

# FIT --------------------------------------------------------------------------

# Compile
compile(here::here("src", "negbin_1re_cumsum.cpp"))
dyn.load(dynlib(here::here("src", "negbin_1re_cumsum")))

ssdr_list <- vector(length = length(fishery_list), mode = "list")
for (i in seq_along(fishery_list)) {
  dum <- fishery_list[[i]]
  obj <- MakeADFun(dum$data, dum$parameters, random = c("z1_k"), 
                  DLL = "negbin_1re_cumsum")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  f_name <- paste(dum$data_type, "negbin_spline_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generated_data", "model_fits", f_name))

  ssdr_list[[i]] <- ssdr
}
  
ssdr_list <- map(fishery_list, function(x) {
  f_name <- paste(x$data_type, "negbin_ssdr.RDS", sep = "_")
  readRDS(here::here("generated_data", "model_fits", f_name))
})

# color pallette
pal <- readRDS(here::here("generated_data", "color_pal.RDS"))  

pred_plot_list <- map2(ssdr_list, fishery_list, function(in_ssdr, in_list) {
  catch <- in_list$input_data
  ssdr <- in_ssdr
  data_type <- in_list$data_type
  
  pred_catch <- in_list$pred_data
  agg_pred_catch <- pred_catch %>% 
    select(-area) %>% 
    distinct()
  
  # plot simple gam predictions for comparison  
  # pred_y <- predict.gam(in_list$m1, 
  #                       pred_catch, type = "response") 
  # pred_catch2 <- pred_catch %>% 
  #   mutate(y = pred_y)
  
  # abundance across areas
  ylab = ifelse(data_type == "rec", "Predicted Monthly Catch Rate", 
                "Predicted Daily Catch Rate")
  abund_pred <- ssdr[rownames(ssdr) %in% "pred_abund", ] 
  pred_ci <- data.frame(pred_est = abund_pred[ , "Estimate"],
                        pred_se =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_catch, .) %>% 
    mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
           pred_up = pred_est + (qnorm(0.975) * pred_se),
           area = fct_reorder(area, as.numeric(region)))
  
  area_preds <- ggplot(data = pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est, colour = region )) +
    geom_ribbon(aes(ymin = pred_low, ymax = pred_up, fill = region), 
                alpha = 0.5)  +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    # geom_line(data = pred_catch2, aes(x = month_n, y = y), col = "black") +
    facet_wrap(~area, scales = "free_y") +
    labs(x = "Month", y = ylab) +
    ggsidekick::theme_sleek()
  
  # cumulative abundance across regions
  agg_abund_pred <- ssdr[rownames(ssdr) %in% "agg_pred_abund", ] 
  agg_pred_ci <- data.frame(pred_est = agg_abund_pred[ , "Estimate"],
                            pred_se =  agg_abund_pred[ , "Std. Error"]) %>%
    cbind(agg_pred_catch, .) %>% 
    mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
           pred_up = pred_est + (qnorm(0.975) * pred_se)) 
  
  agg_preds <- ggplot(data = agg_pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est, colour = region)) +
    geom_ribbon(aes(ymin = pred_low, ymax = pred_up, fill = region), 
                alpha = 0.5) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    labs(x = "Month", y = ylab) +
    ggsidekick::theme_sleek()
  
  f_name <- paste(data_type, "negbin_spline_prediction.pdf", sep = "_")
  pdf(here::here("figs", "model_pred", "neg_bin_only", f_name))
  print(area_preds)
  print(agg_preds)
  dev.off()
  
  f_name2 <- paste(data_type, "negbin_spline_predictions.RDS", sep = "_")
  plot_list_out <- list(area_preds, agg_preds)
  saveRDS(plot_list_out, 
          here::here("figs", "model_pred", "neg_bin_only", f_name2))
})
