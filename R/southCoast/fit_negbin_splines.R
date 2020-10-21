## Neg bin model fit w/ splines
# July 20, 2020
# Fit neg binomial model to catch data from WCVI troll fishery and southern BC
# rec fisheries
# Modified from fit_negbin.R to include splines

library(tidyverse)
library(TMB)
library(mgcv)

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

# ggplot(rec_catch, aes(x = log(eff), y = catch, fill = region)) +
#   geom_point(shape = 21) +
#   # stat_smooth(method = "gam", k = 2) +
#   facet_wrap(~area_n)
# ggplot(comm_catch, aes(x = log(eff), y = catch, fill = region)) +
#   geom_point(shape = 21) +
#   # stat_smooth(method = "gam", k = 2) +
#   facet_wrap(~area_n)
# ggplot(rec_catch) +
#   geom_boxplot(aes(x = area, y = catch, fill = region)) 
# ggplot(rec_catch) +
#   geom_boxplot(aes(x = area, y = cpue, fill = region)) 
# 
# rec_catch %>% 
#   group_by(region, month, year) %>% 
#   summarize(sum_catch = sum(catch)) %>% 
#   ggplot(.) +
#   geom_boxplot(aes(x = month, y = sum_catch, fill = region)) +
#   facet_wrap(~region)
# comm_catch %>% 
#   group_by(region, month, year) %>% 
#   summarize(sum_catch = sum(catch)) %>% 
#   ggplot(.) +
#   geom_boxplot(aes(x = month, y = sum_catch, fill = region)) +
#   facet_wrap(~region)
# 
# ggplot(rec_catch) +
#   geom_histogram(aes(x = eff))
# ggplot(comm_catch) +
#   geom_histogram(aes(x = eff))

# PREP DATA --------------------------------------------------------------------

catch_dat <- rec_catch

prep_catch <- function (catch_dat, data_type = NULL) {
  yr_vec <- as.numeric(as.factor(as.character(catch_dat$year))) - 1
  
  #generate model matrix based on GAM
  months <- unique(catch_dat$month_n)
  n_months <- length(months)
  spline_type <- ifelse(n_months == 12, "cc", "tp")
  n_knots <- ifelse(max(months) == 12, 4, 3)
  
  m1 <- gam(catch ~ area + s(month_n, bs = spline_type, k = n_knots, by = area) +
              offset,
            knots = list(month_n = c(min(months), max(months))),
            data = catch_dat,
            family = nb)
  fix_mm <- predict(m1, type = "lpmatrix")
  offset_pos <- grep("^offset$", colnames(fix_mm))
  
  # make predictive model matrix including null values for effort
  pred_dat <- expand.grid(
    month_n = seq(min(catch_dat$month_n), 
                  max(catch_dat$month_n),
                  length.out = 50),
    area = unique(catch_dat$area),
    offset = mean(catch_dat$offset)
  ) %>% 
    left_join(., catch_dat %>% select(region, region_c, area), by = "area") %>% 
    mutate(month = as.factor(round(month_n, 3)),
           order = row_number(),
           reg_month = paste(region, month, sep = "_"),
           reg_month_f = fct_reorder(factor(reg_month), order)) %>% 
    select(-order, -reg_month) %>% 
    distinct()
  pred_mm_catch <- predict(m1, pred_dat, type = "lpmatrix")
  
  # construct factor key for regional aggregates associated with areas
  grouping_vec <- as.numeric(pred_dat$reg_month_f) - 1
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
  
  # fix starting value of offset to one, then map
  parameters$b1_j[offset_pos] <- 1
  b_j_map <- seq_along(parameters$b1_j)
  b_j_map[offset_pos] <- NA
  tmb_map <- list(b1_j = as.factor(b_j_map))
  
  if (is.null(data_type)) {
    data_type <- unique(catch_dat$legal)
  }
  
  list("fix_mm" = fix_mm, "mm_pred" = pred_mm_catch, "pred_data" = pred_dat, 
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "input_data" = catch_dat, "m1" = m1, "tmb_map" = tmb_map)
}

comm_list <- prep_catch(comm_catch, data_type = "comm")
rec_list <- prep_catch(rec_catch, data_type = "rec")
fishery_list <- list(comm_list, rec_list)

# FIT --------------------------------------------------------------------------

# Compile
compile(here::here("src", "negbin_1re_cumsum_rw.cpp"))
dyn.load(dynlib(here::here("src", "negbin_1re_cumsum_rw")))

for (i in seq_along(fishery_list)) {
  dum <- fishery_list[[i]]
  obj <- MakeADFun(data = dum$data, parameters = dum$parameters, 
                   random = c("z1_k"), 
                   DLL = "negbin_1re_cumsum_rw", map = dum$tmb_map)
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  f_name <- paste(dum$data_type, "negbin_spline_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generated_data", "model_fits", f_name))
}
  
ssdr_list <- map(fishery_list, function(x) {
  f_name <- paste(x$data_type, "negbin_spline_ssdr.RDS", sep = "_")
  readRDS(here::here("generated_data", "model_fits", f_name))
})

# color pallette
pal <- readRDS(here::here("generated_data", "color_pal.RDS"))  

pred_plot_list <- map2(ssdr_list, fishery_list, function(in_ssdr, in_list) {
  # in_list <- fishery_list[[2]]
  # ssdr <- ssdr_list[[2]]
  ssdr <- in_ssdr
  catch <- in_list$input_data
  data_type <- in_list$data_type
  pred_catch <- in_list$pred_data
  
  # Check data fit
  fit_dat <- in_list$input_data %>% 
    select(region_c, area, month_n, year, catch, eff, offset) %>% 
    mutate(log_catch = log(catch) - offset,
           log_cpue = log(catch) / offset,
           obs_cpue = exp(log_cpue),
           mu_log_catch = as.numeric(ssdr[rownames(ssdr) %in% "s1", "Estimate"]),
           mu_log_catch2 = mu_log_catch - offset,
           mu_log_cpue = mu_log_catch / mean(offset))
  fit_plot <- ggplot(fit_dat) +
    geom_point(aes(x = month_n, y = mu_log_catch2), colour = "black") +
    geom_point(aes(x = month_n, y = log_catch), colour = "red") +
    facet_grid(area~year, scales = "free_y") +
    ggsidekick::theme_sleek()
  
  # abundance across areas
  ylab = "Predicted Monthly Catch Rate"
  abund_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] 
  pred_ci <- data.frame(est_link = abund_pred[ , "Estimate"],
                        se_link =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_catch, .) %>% 
    mutate(
      pred_est = exp(est_link),
      pred_low = exp(est_link + (qnorm(0.025) * se_link)),
      pred_up = exp(est_link + (qnorm(0.975) * se_link)),
      pred_est_logcpue = est_link / offset,
      pred_low_logcpue = (est_link + (qnorm(0.025) * se_link)) / offset,
      pred_up_logcpue = (est_link + (qnorm(0.975) * se_link)) / offset
    )

  log_area_preds <- ggplot(data = pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est_logcpue, colour = region_c)) +
    geom_ribbon(aes(ymin = pred_low_logcpue, ymax = pred_up_logcpue, 
                    fill = region_c), 
                alpha = 0.5) +
    geom_point(data = fit_dat, aes(x = month_n, y = log_cpue, fill = region_c),
               shape = 21) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    scale_x_continuous(breaks = seq(1, 12, by = 1), limits = c(1, 12)) +
    facet_wrap(~area, scales = "free_y") +
    labs(x = "Month", y = ylab) +
    ggsidekick::theme_sleek()
  
  # random intercepts
  rand_int <- ssdr[rownames(ssdr) %in% "z1_k", ] 
  year_ints <- data.frame(year = levels(as.factor(as.character(catch$year))),
                          z1_k = as.numeric(rand_int[, "Estimate"])) 
  year_preds <- expand.grid(year = year_ints$year, 
                            month = unique(pred_catch$month)) %>% 
    left_join(pred_catch %>% select(area, month), ., by = "month") %>% 
    left_join(., year_ints, by = "year")
  
  yr_pred_ci <- pred_ci %>% 
    select(month_n:est_link) %>% 
    left_join(., year_preds, by = c("area", "month")) %>% 
    mutate(
      year = as.factor(year),
      est_link_yr = est_link + z1_k,
      pred_est_logcpue = est_link_yr / offset
    ) 
  
  yr_log_area_preds <- ggplot(data = yr_pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est_logcpue, colour = year)) +
    geom_point(data = fit_dat, aes(x = month_n, y = log_cpue, fill = year),
               shape = 21) +
    scale_fill_discrete(name = "Year") +
    scale_colour_discrete(name = "Year") +
    scale_x_continuous(breaks = seq(1, 12, by = 1), limits = c(1, 12)) +
    facet_wrap(~area, scales = "free_y") +
    labs(x = "Month", y = ylab) +
    ggsidekick::theme_sleek()
    
  # cumulative abundance across regions
  agg_pred_catch <- pred_catch %>% 
    select(-area) %>% 
    distinct()
  n_reg <- fit_dat %>% 
    group_by(region_c) %>% 
    summarize(n_areas = length(unique(area)),
              .groups = "drop") %>% 
    ungroup()
  agg_abund_pred <- ssdr[rownames(ssdr) %in% "log_agg_pred_abund", ] 
  agg_pred_ci <- data.frame(est_link = agg_abund_pred[ , "Estimate"],
                            se_link =  agg_abund_pred[ , "Std. Error"]) %>%
    cbind(agg_pred_catch, .) %>% 
    left_join(., n_reg, by = "region_c") %>% 
    mutate(
      pred_est = exp(est_link),
      pred_low = exp(est_link + (qnorm(0.025) * se_link)),
      pred_up = exp(est_link + (qnorm(0.975) * se_link)),
      pred_est_logcpue = est_link / log(n_areas*exp(offset)),
      pred_low_logcpue = (est_link + (qnorm(0.025) * se_link)) / log(n_areas*exp(offset)),
      pred_up_logcpue = (est_link + (qnorm(0.975) * se_link)) / log(n_areas*exp(offset))
    ) 
  
  agg_obs_dat <- fit_dat %>% 
    group_by(region_c, month_n, year) %>% 
    summarize(sum_catch = sum(catch),
              sum_eff = sum(eff),
              .groups = "drop") %>% 
    ungroup() %>% 
    mutate(log_agg_cpue = log(sum_catch) / log(mean(sum_eff)))
  
  log_agg_preds <- ggplot(data = agg_pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est_logcpue, colour = region_c)) +
    geom_ribbon(aes(ymin = pred_low_logcpue, ymax = pred_up_logcpue, 
                    fill = region_c), 
                alpha = 0.5) +
    geom_point(data = agg_obs_dat, aes(x = month_n, y = log_agg_cpue,
                                       fill = region_c),
               shape = 21) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_x_continuous(breaks = seq(1, 12, by = 1), limits = c(1, 12)) +
    scale_colour_manual(name = "Region", values = pal) +
    labs(x = "Month", y = ylab) +
    ggsidekick::theme_sleek()
  
  # aggregate year-specific effects
  agg_year_preds <- expand.grid(year = year_ints$year, 
                            month = unique(pred_catch$month)) %>% 
    left_join(agg_pred_catch %>% select(region, month), ., 
              by = "month") %>% 
    left_join(., year_ints, by = "year")
  
  yr_agg_pred_ci <- agg_pred_ci %>% 
    select(month_n:est_link) %>% 
    left_join(., agg_year_preds, by = c("region", "month")) %>% 
    mutate(
      year = as.factor(year),
      est_link_yr = est_link + z1_k,
      pred_est_logcpue = est_link_yr / offset
    ) 
  
  yr_log_agg_preds <- ggplot(data = yr_agg_pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est_logcpue, colour = year)) +
    geom_point(data = agg_fit_dat, aes(x = month_n, y = log_agg_cpue2, fill = year),
               shape = 21) +
    scale_fill_discrete(name = "Year") +
    scale_colour_discrete(name = "Year") +
    scale_x_continuous(breaks = seq(1, 12, by = 1), limits = c(1, 12)) +
    facet_wrap(~region_c, scales = "free_y") +
    labs(x = "Month", y = ylab) +
    ggsidekick::theme_sleek()
  
  f_name <- paste(data_type, "negbin_spline_prediction.pdf", sep = "_")
  pdf(here::here("figs", "model_pred", "neg_bin_only", f_name))
  print(area_preds)
  print(agg_preds)
  dev.off()
  
  f_name2 <- paste(data_type, "negbin_spline_predictions.RDS", sep = "_")
  plot_list_out <- list(area_preds, agg_preds, pred_ci, agg_pred_ci)
  saveRDS(plot_list_out, 
          here::here("figs", "model_pred", "neg_bin_only", f_name2))
})

