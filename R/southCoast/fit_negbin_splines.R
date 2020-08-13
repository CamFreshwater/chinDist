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
    mutate(region_c = as.character(region),
           region = factor(abbreviate(region, minlength = 4)),
           area_n = as.numeric(area),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff)),
           offset = log(eff)
    ) %>% 
    arrange(region, month) 
}

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  # calculate monthly catch and effort 
  group_by(catchReg, area, month, year) %>% 
  summarize(catch = sum(catch), 
            eff = sum(boatDays), 
            cpue = catch / eff,
            .groups = "drop") %>% 
  rename(region = catchReg) %>% 
  clean_catch(.) %>% 
  # drop inside areas where seasonal catches not available
  filter(!area_n < 100) %>% 
  droplevels() %>% 
  select(region, region_c, area, area_n, month, month_n, year, catch, eff, cpue, 
         eff_z, offset)

# tally_f <- function(dat) {
#   dat %>% 
#     group_by(region, area) %>% 
#     tally()
# }
# tally_f(comm_catch) %>% 
#   print(n = Inf)
# tally_f(comm_catch2) %>% 
#   print(n = Inf)

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
            cpue = catch / eff,
            .groups = "drop") %>% 
  clean_catch(.) %>% 
  # drop months with minimal catch estimates
  filter(!month_n < 5,
         !month_n > 9) %>% 
  # drop areas with fewer than 10 (i.e annual estimate per subarea = 1 sample)
  group_by(area_n) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(!n < 10) %>% 
  droplevels()

table(rec_catch$month_n, rec_catch$area, rec_catch$region)


ggplot(rec_catch, aes(x = log(eff), y = catch, fill = region)) +
  geom_point(shape = 21) +
  # stat_smooth(method = "gam", k = 2) +
  facet_wrap(~area_n)
ggplot(comm_catch, aes(x = log(eff), y = catch, fill = region)) +
  geom_point(shape = 21) +
  # stat_smooth(method = "gam", k = 2) +
  facet_wrap(~area_n)
ggplot(rec_catch) +
  geom_boxplot(aes(x = area, y = catch, fill = region)) 
ggplot(rec_catch) +
  geom_boxplot(aes(x = area, y = cpue, fill = region)) 

rec_catch %>% 
  group_by(region, month, year) %>% 
  summarize(sum_catch = sum(catch)) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = month, y = sum_catch, fill = region)) +
  facet_wrap(~region)
comm_catch %>% 
  group_by(region, month, year) %>% 
  summarize(sum_catch = sum(catch)) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = month, y = sum_catch, fill = region)) +
  facet_wrap(~region)

ggplot(rec_catch) +
  geom_histogram(aes(x = eff))
ggplot(comm_catch) +
  geom_histogram(aes(x = eff))

# PREP DATA --------------------------------------------------------------------

catch_dat <- comm_catch

prep_catch <- function (catch_dat, data_type = NULL, n_knots = 4) {
  yr_vec <- as.numeric(as.factor(as.character(catch_dat$year))) - 1
  
  #generate model matrix based on GAM
  months <- unique(catch_dat$month_n)
  n_months <- length(months)
  spline_type <- ifelse(n_months == 12, "cc", "tp")
  # m1 <- gam(catch ~ area + s(month_n, bs = spline_type, k = n_knots, by = area) +
  #              s(eff_z, bs = "tp", k = 4),
  #            knots = list(month_n = c(min(months), max(months))),
  #            data = catch_dat,
  #            family = nb)
  m1 <- gam(catch ~ area + s(month_n, bs = spline_type, k = n_knots, by = area) +
              offset,
            knots = list(month_n = c(min(months), max(months))),
            data = catch_dat,
            family = nb)
  fix_mm <- predict(m1, type = "lpmatrix")
  # fix_mm2 <- predict(m2, type = "lpmatrix")
  offset_pos <- grep("^offset$", colnames(fix_mm))
  
  # make predictive model matrix including null values for effort
  pred_dat <- expand.grid(
    month_n = seq(min(catch_dat$month_n), 
                  max(catch_dat$month_n),
                  length.out = 50),
    area = unique(catch_dat$area),
    offset = mean(catch_dat$offset)
    # eff_z = 0
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
compile(here::here("src", "negbin_1re_cumsum.cpp"))
dyn.load(dynlib(here::here("src", "negbin_1re_cumsum")))

for (i in seq_along(fishery_list)) {
  dum <- fishery_list[[i]]
  obj <- MakeADFun(dum$data, dum$parameters, random = c("z1_k"), 
                  DLL = "negbin_1re_cumsum", map = dum$tmb_map)
  
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
  # in_list <- fishery_list[[1]]
  # ssdr <- ssdr_list[[2]]
  ssdr <- in_ssdr
  catch <- in_list$input_data
  data_type <- in_list$data_type
  
  pred_catch <- in_list$pred_data
  agg_pred_catch <- pred_catch %>% 
    select(-area) %>% 
    distinct()
  
  # abundance across areas
  ylab = ifelse(data_type == "rec", "Predicted Monthly Catch Rate", 
                "Predicted Daily Catch Rate")
  abund_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] 
  pred_ci <- data.frame(est_link = abund_pred[ , "Estimate"],
                        se_link =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_catch, .) %>% 
    mutate(
      pred_est = exp(est_link),
      pred_low = exp(est_link + (qnorm(0.025) * se_link)),
      pred_up = exp(est_link + (qnorm(0.975) * se_link))
    ) 

  area_preds <- ggplot(data = pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est, colour = region_c)) +
    geom_ribbon(aes(ymin = pred_low, ymax = pred_up, fill = region_c), 
                alpha = 0.5)  +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    facet_wrap(~area, scales = "free_y") +
    labs(x = "Month", y = ylab) +
    ggsidekick::theme_sleek()
  
  # cumulative abundance across regions
  agg_abund_pred <- ssdr[rownames(ssdr) %in% "log_agg_pred_abund", ] 
  agg_pred_ci <- data.frame(est_link = agg_abund_pred[ , "Estimate"],
                            se_link =  agg_abund_pred[ , "Std. Error"]) %>%
    cbind(agg_pred_catch, .) %>% 
    mutate(
      pred_est = exp(est_link),
      pred_low = exp(est_link + (qnorm(0.025) * se_link)),
      pred_up = exp(est_link + (qnorm(0.975) * se_link))
    ) 
  
  agg_preds <- ggplot(data = agg_pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_est, colour = region_c)) +
    geom_ribbon(aes(ymin = pred_low, ymax = pred_up, fill = region_c), 
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
  plot_list_out <- list(area_preds, agg_preds, pred_ci, agg_pred_ci)
  saveRDS(plot_list_out, 
          here::here("figs", "model_pred", "neg_bin_only", f_name2))
})

