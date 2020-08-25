## Dirichlet model fit
# July 9, 2020
# Fit dirichlet model to stock composition data at various aggregate levels
# Accounts for uncertain GSI and uses julian day within a strata as a sampling 
# event

library(tidyverse)
library(TMB)
library(ggplot2)

# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(legal == "legal") %>% 
  mutate(sample_id = paste(temp_strata, jDay, year, sep = "_"))
# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds")) %>%
  mutate(sample_id = paste(temp_strata, jDay, year, sep = "_"))
  

# helper function to calculate aggregate probs
calc_agg_prob <- function(grouped_data, full_data) {
  grouped_data %>% 
    summarize(agg_prob = sum(adj_prob)) %>% 
    arrange(sample_id, desc(agg_prob)) %>%
    ungroup() %>% 
    distinct() %>% 
    left_join(.,
              full_data %>% 
                select(sample_id, temp_strata:area_n) %>% 
                distinct(),
              by = "sample_id")
}

# helper function to pool non-focal stocks
pool_aggs <- function(full_data) {
  full_data %>% 
    mutate(
      pst_agg = case_when(
        grepl("CR-", pst_agg) ~ "CR",
        grepl("CST", pst_agg) ~ "CA/OR/WA",
        pst_agg %in% c("NBC_SEAK", "WCVI") ~ "BC-coast",
        pst_agg %in% c("PSD", "SOG") ~ "SalSea",
        TRUE ~ pst_agg
      ),
      reg1 = case_when(
        Region1Name %in% c("Fraser_Fall", "ECVI", "Fraser_Summer_4.1", 
                           "Fraser_Spring_4.2", "Fraser_Spring_5.2", 
                           "Fraser_Summer_5.2", "WCVI", "SOMN") ~ Region1Name,
        TRUE ~ "Other"
      )
    )
}

# helper function to use above for gsi and cwt samples
clean_comp <- function(grouping_col, raw_data, ...) {
  raw_data %>% 
    pool_aggs() %>% 
    rename(agg = grouping_col) %>%
    group_by(sample_id, agg) %>%
    calc_agg_prob(., raw_data)
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
    data = pmap(list(grouping_col, raw_data), .f = clean_comp)
  )

## PREP INPUTS ----------------------------------------------------------------

comp_in <- comp$data[[3]]

prep_dir_inputs <- function(comp_in, data_type) {
  gsi_trim <- comp_in %>% 
    group_by(region, month, year) %>%
    mutate(nn = sum(agg_prob)) %>%
    #remove strata with less than 100 individuals total (necessary for 
    #convergence for Canadian specific predictions)
    filter(!nn < 100) %>%
    ungroup() %>%
    droplevels() %>%
    select(sample_id, region, year, month, season, agg, agg_prob) %>%
    distinct()
  
  gsi_wide <- gsi_trim %>% 
    pivot_wider(., names_from = agg, values_from = agg_prob) %>%
    mutate_if(is.numeric, ~replace_na(., 0.000001))
  
  y_obs <- gsi_wide %>% 
    select(-c(sample_id:season)) %>% 
    as.matrix() 

  yr_vec <- as.numeric(gsi_wide$year) - 1
  fix_mm <- model.matrix(~ region + month, gsi_wide) #fixed covariates only
  
  # data frame for predictions
  pred_dat <- gsi_wide %>%
    select(region, month) %>%
    distinct() %>%
    arrange(region, month)
  pred_mm <- model.matrix(~ region + month, pred_dat)
  
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
  mutate(model_inputs = map2(data, dataset, .f = prep_dir_inputs))
         
# FIT --------------------------------------------------------------------------

compile(here::here("src", "dirichlet_randInt.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_randInt")))

fit_model <- function(x) {
  ## Make a function object
  # x <- comp2$model_inputs[[3]]
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
  
  f_name <- paste(x$data_type, "dirichlet_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generated_data", "model_fits", f_name))
  
  return(ssdr)
}

# fit
map(comp2$model_inputs, .f = fit_model)


## PLOT PREDICTIONS ------------------------------------------------------------

ssdr_list <- map(comp2$model_inputs, function (x) {
  f_name <- paste(x$data_type, "dirichlet_ssdr.RDS", sep = "_")
  readRDS(here::here("generated_data", "model_fits", f_name))
})

pal <- readRDS(here::here("generated_data", "color_pal.RDS"))

plot_list <- map2(comp2$model_inputs, ssdr_list, function(x, ssdr) {
  y_obs <- x$data$y_obs
  k <- ncol(y_obs) # number of stocks
  stk_names <- colnames(y_obs)
  N <- nrow(y_obs)
  pred_dat <- x$pred_dat
  
  pred_dat2 <- purrr::map_dfr(seq_len(k), ~pred_dat)
  pred_mat <- ssdr[rownames(ssdr) %in% "pred_pi_prop", ]
  pred_ci <- data.frame(stock = as.character(rep(stk_names, 
                                                 each = nrow(pred_dat))),
                        pred_prob_est = pred_mat[ , "Estimate"],
                        pred_prob_se =  pred_mat[ , "Std. Error"]) %>% 
    cbind(pred_dat2, .) %>%
    mutate(pred_prob_low = pred_prob_est + (qnorm(0.025) * pred_prob_se),
           pred_prob_up = pred_prob_est + (qnorm(0.975) * pred_prob_se), 
           region = abbreviate(region, minlength = 4)) 
  
  # calculate raw proportion data for comparison
  raw_prop <- x$long_data %>% 
    group_by(region, month, year, agg) %>%
    summarize(samp_g = sum(agg_prob)) %>% 
    group_by(region, month, year) %>%
    mutate(samp_total = sum(samp_g)) %>% 
    ungroup() %>% 
    mutate(samp_g_ppn = samp_g / samp_total,
           stock = fct_reorder(agg, desc(samp_g_ppn)), 
           region = abbreviate(region, minlength = 4)) 
  
  pred_plot <- ggplot() +
    geom_pointrange(data = pred_ci,
                    aes(x = as.numeric(as.character(month)), y = pred_prob_est, 
                        ymin = pred_prob_low, ymax = pred_prob_up, 
                        fill = region), 
                    shape = 21, size = 0.4, position = position_dodge(0.6)) +
    labs(y = "Probability", x = "Month") +
    scale_x_continuous(breaks = seq(1, 12, by = 1)) +
    scale_fill_manual(name = "Region", values = pal) +
    facet_wrap(~stock) +
    ggsidekick::theme_sleek()
  pred_plot_raw <- pred_plot +
    geom_point(data = raw_prop, 
               aes(x = as.numeric(as.character(month)), 
                   y = samp_g_ppn, fill = region),
               shape = 21, alpha = 0.4, position = position_dodge(0.6)) 
  
  f_name <- paste(x$data_type, "dirichlet_pred.pdf", sep = "_")
  pdf(here::here("figs", "model_pred", "dirichlet_only", f_name))
  print(pred_plot_raw)
  dev.off()
  
  f_name2 <- paste(x$data_type, "dirichlet_pred.rds", sep = "_")
  saveRDS(pred_plot,
          here::here("figs", "model_pred", "dirichlet_only", f_name2))
})
