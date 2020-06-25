## Neg bin model fit
# June 2, 2020
# Fit neg binomial model to catch data from WCVI troll fishery and southern BC
# rec fisheries
# Changed June 11 so that fit to areas, then predictions summed

library(tidyverse)
library(TMB)

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
  # drop month where catch data are missing even though gsi samples available
  mutate(temp_strata = paste(month_n, region, sep = "_")) %>% 
  filter(!temp_strata == "7_SWVI") 

#recreational catch data
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "monthlyCatch_rec.RDS")) %>% 
  #group by subarea to get rid of adipose and released legal categories
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
  clean_catch(.) #%>% 
  # filter(month %in% c("6", "7", "8")) %>% 
  # droplevels()

table(rec_catch$area, rec_catch$month)
table(comm_catch$area, comm_catch$month)


# visualize relative catch/effort
temp <- comm_catch %>% 
  select(month, year, region, eff, catch, cpue) %>% 
  mutate(dat = "comm")
plot_catch <- rec_catch %>% 
  mutate(dat = legal) %>% 
  select(month, year, region, eff, catch, cpue, dat) %>% 
  rbind(., temp)

ggplot(plot_catch) +
  geom_histogram(aes(x = cpue)) +
  facet_wrap(~dat, scales = "free")
ggplot(plot_catch, aes(x = eff, y = catch)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam") +
  facet_wrap(~dat, scales = "free") +
  ggsidekick::theme_sleek()

rec_catch %>% 
  filter(legal == "legal") %>% 
  ggplot(.) + 
  geom_boxplot(aes(x = area, y = cpue)) +
  facet_wrap(~region, scales = "free_x")


# PREP DATA --------------------------------------------------------------------

# function to make predictive dataframes 
make_pred_dat <- function (dat) {
  dat %>%
    select(area, month, region) %>%
    distinct() %>%
    arrange(area, month, region)
}

prep_catch <- function (catch_dat, effort ="overall", data_type = NULL) {
  yr_vec <- as.numeric(as.factor(as.character(catch_dat$year))) - 1
  
  # model matrix for fixed effects
  fix_mm <- model.matrix(~ area + month + eff_z + eff_z2, catch_dat)
  
  # use average effort by month and area for predictions
  if (effort == "month") {
    mean_eff <- catch_dat %>% 
      group_by(area, month) %>%
      summarize(eff_z = mean(eff_z),
                eff_z2 = mean(eff_z2)) %>% 
      ungroup()
    pred_dat_catch <- make_pred_dat(catch_dat) %>% 
      left_join(., mean_eff, by = c("area", "month")) %>% 
      mutate(reg_month = paste(region, month, sep = "_"))
    pred_mm_catch <- model.matrix(~ area + month + eff_z + eff_z2, 
                                  pred_dat_catch) 
  }  
  # set effort to zero for all predictions
  if (effort == "overall") {
    pred_dat_catch <- make_pred_dat(catch_dat) %>% 
      mutate(reg_month = paste(region, month, sep = "_"))
    pred_mm_catch <- model.matrix(~ area + month, pred_dat_catch) %>%
      cbind(.,
            eff_z = rep(0, n = nrow(.)),
            eff_z2 = rep(0, n = nrow(.)))
  }
  
  # construct factor key for regional aggregates associated with areas
  grouping_vec <- as.numeric(as.factor(as.character(pred_dat_catch$reg_month))) - 1
  grouping_key <- unique(grouping_vec)
  
  
  data <- list(y1_i = catch_dat$catch,
               X1_ij = fix_mm,
               factor1k_i = yr_vec,
               nk1 = length(unique(yr_vec)),
               X1_pred_ij = pred_mm_catch,
               pred_factor2k_h = grouping_vec,
               pred_factor2k_levels = grouping_key
  )
  
  # Fit simple model to initialize tmb 
  m1 <- lm(log(catch + 0.0001) ~ area + month + eff_z + eff_z2, 
           data = catch_dat)
  
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
  
  list("fix_mm" = fix_mm, "mm_pred" = pred_mm_catch, 
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "input_data" = catch_dat)
}

comm_list <- prep_catch(comm_catch, data_type = "comm")

# for recreational data prep separate inputs 
rec_list1 <- rec_catch %>% 
  split(., .$legal) 
nms <- names(rec_list1)
rec_list <- rec_list1 %>%
  map(., prep_catch)

fishery_list <- list(comm_list, rec_list[[1]], rec_list[[2]])
names(fishery_list) <- c("comm", nms)


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
  
  f_name <- paste(dum$data_type, "negbin_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generated_data", "model_fits", f_name))
  
  ssdr_list[[i]] <- ssdr
}


# PREDICTIONS ------------------------------------------------------------------

ssdr_list <- map(fishery_list, function(x) {
  f_name <- paste(x$data_type, "negbin_ssdr.RDS", sep = "_")
  readRDS(here::here("generated_data", "model_fits", f_name))
})

pal <- readRDS(here::here("generated_data", "color_pal.RDS"))

pred_plot_list <- map2(ssdr_list, fishery_list, function(in_ssdr, in_list) {
  catch <- in_list$input_data
  ssdr <- in_ssdr
  # fac_key <- in_list$fac_key
  data_type <- in_list$data_type

  pred_catch <- make_pred_dat(catch) %>% 
    mutate(reg_month = paste(region, month, sep = "_"))
  
  abund_pred <- ssdr[rownames(ssdr) %in% "pred_abund", ] 
  pred_ci <- data.frame(pred_est = abund_pred[ , "Estimate"],
                        pred_se =  abund_pred[ , "Std. Error"]) %>%
    cbind(pred_catch, .) %>% 
    mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
           pred_up = pred_est + (qnorm(0.975) * pred_se)) 
  real_preds <- ggplot() +
    geom_pointrange(data = pred_ci, 
                    aes(x = as.factor(month), y = pred_est, ymin = pred_low, 
                        ymax = pred_up, fill = region), shape = 21) +
    scale_fill_manual(name = "Region", values = pal) +
    labs(x = "Month", y = "Predicted Real Catch (mean effort)",
         title = data_type) +
    ggsidekick::theme_sleek() +
    facet_wrap(~ fct_reorder(area, as.numeric(region)))
  
  
  agg_abund_pred <- ssdr[rownames(ssdr) %in% "agg_pred_abund", ] 
  agg_pred_ci <- data.frame(pred_est = agg_abund_pred[ , "Estimate"],
                            pred_se =  agg_abund_pred[ , "Std. Error"]) %>%
    cbind(pred_catch %>% 
            select(region, month, reg_month) %>% 
            distinct(), 
          .) %>% 
    mutate(pred_low = pred_est + (qnorm(0.025) * pred_se),
           pred_up = pred_est + (qnorm(0.975) * pred_se))
  agg_preds <- ggplot() +
    geom_pointrange(data = agg_pred_ci, 
                    aes(x = as.factor(month), y = pred_est,
                        ymin = pred_low, ymax = pred_up, fill = region), 
                    shape = 21) +
    scale_fill_manual(name = "Region", values = pal) +
    labs(x = "month", y = "predicted real catch (mean effort)",
         title = data_type) +
    ggsidekick::theme_sleek() +
    facet_wrap(~region)
  
  ## plot predictions across different levels of effort
  # pull coeficients
  abund_b <- ssdr[rownames(ssdr) %in% "b1_j", ]
  
  # generate data for each factor level
  pred_eff <- catch %>% 
    select(area, month, eff_z, eff_z2) %>% 
    group_by(area, month) %>% 
    sample_n(., size = 50, replace = TRUE) %>% 
    ungroup()
  pred_mm2 <- model.matrix(~ area + month + eff_z + eff_z2, pred_eff)
  pred_catch <- pred_mm2 %*% abund_b 
  pred_plot <- pred_eff %>% 
    mutate(log_catch = pred_catch[ , 'Estimate'],
           catch = exp(log_catch),
           dataset = "pred")
  
  if (length(unique(catch$area)) > 9) {
    subset_areas <- sample(unique(catch$area), size = 9)
    subset_months <- sample(unique(catch$month), size = 3)
  } else {
    subset_areas <- unique(catch$area)
    subset_months <- unique(catch$month)
  }
  
  var_effort_preds <- catch %>% 
    mutate(dataset = "obs") %>% 
    select(catch, dataset, area, month) %>% 
    rbind(., 
          pred_plot %>% 
            select(catch, dataset, area, month)
    ) %>% 
    filter(area %in% subset_areas,
           month %in% subset_months) %>% 
    ggplot(.) +
    geom_boxplot(aes(x = month, y = catch, fill = dataset)) +
    facet_wrap(~ area, nrow = 3, scales = "free_y") +
    ggsidekick::theme_sleek() +
    labs(fill = "Data", title = data_type)
  
  f_name <- paste(data_type, "negbin_prediction.pdf", sep = "_")
  pdf(here::here("figs", "model_pred", "neg_bin_only", f_name))
  print(real_preds)
  print(agg_preds)
  print(var_effort_preds)
  dev.off()
  
  f_name2 <- paste(data_type, "negbin_predictions.RDS", sep = "_")
  plot_list_out <- list(real_preds, agg_preds)
  saveRDS(plot_list_out, 
          here::here("figs", "model_pred", "neg_bin_only", f_name2))
})

