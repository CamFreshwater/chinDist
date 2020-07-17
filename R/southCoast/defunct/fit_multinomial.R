## Multinomial model fit
# January 20, 2020
# Fit multinomial model to stock composition data at various aggregate levels
# (assumes perfect GSI assignment); in both cases aggregate probability reflects
# summed probabilities of a given region of origin (ie reg1 or 3) for a given
# individual
# Also plot fixed effect version to evaluate impact on prediction intervals

library(tidyverse)
library(TMB)
library(ggplot2)

# cwt data 
cwt <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                  "cwt_recovery_clean.rds")) %>% 
  select(gear, subset_data) %>% 
  unnest(., cols = c(subset_data)) %>% 
  select(id:date, gear, pres:pst_agg) %>% 
  ungroup()
# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds"))
# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                          "wcviIndProbsLong.rds"))

# helper function to calculate aggregate probs
calc_max_prob <- function(grouped_data, full_data, thresh = 0.75) {
  out <- grouped_data %>% 
    summarize(agg_prob = sum(adj_prob)) %>% 
    arrange(id, desc(agg_prob)) %>%
    ungroup() %>% 
    group_by(id) %>% 
    mutate(max_assignment = max(agg_prob)) %>% 
    # Remove samples where top stock ID is less than threshold probability
    filter(!agg_prob < max_assignment, 
           !max_assignment < thresh) %>% 
    ungroup() %>% 
    distinct() %>% 
    left_join(.,
              full_data %>% 
                select(id:area_n) %>% 
                distinct(),
              by = "id") %>% 
    select(-max_assignment)
  
  # colnames(out)[2] <- "agg"
  
  return(out)
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
clean_comp <- function(sample, grouping_col, raw_data, ...) {
  if (sample == "gsi") {
    out <- raw_data %>% 
      pool_aggs() %>% 
      rename(agg = grouping_col) %>%
      group_by(id, agg) %>%
      calc_max_prob(., raw_data, thresh = 0.75)
  } 
  if (sample == "cwt") {
    out <- raw_data %>%
      pool_aggs() %>% 
      rename(agg = grouping_col) %>% 
      select(id, agg, temp_strata:area_n)
  }
  return(out)
}

#large scale PST groups
# rec_pst <- rec %>% 
#   pool_aggs() %>% 
#   group_by(id, pst_agg) %>%
#   calc_max_prob(., rec, thresh = 0.75)
# comm_pst <- comm %>% 
#   pool_aggs() %>% 
#   group_by(id, pst_agg) %>%
#   calc_max_prob(., comm, thresh = 0.75)
# #aggregate with Canadian CU focus
# rec_can <- rec %>% 
#   pool_aggs() %>% 
#   group_by(id, reg1) %>%
#   calc_max_prob(., rec, thresh = 0.75)
# comm_can <- comm %>% 
#   pool_aggs() %>% 
#   group_by(id, reg1) %>%
#   calc_max_prob(., comm, thresh = 0.75)
# 
comm_cwt <- cwt %>% filter(gear == "troll")
rec_cwt <- cwt %>% filter(gear == "sport")

# combined tibble 
comp <- tibble(
  sample = c(rep("gsi", 4), rep("cwt", 4)),
  fishery = rep(c("troll", "sport"), times = 4),
  grouping = rep(c(rep("pst", 2), rep("can", 2)), 2),
  dataset = paste(sample, fishery, grouping, sep = "_"),
  # incorporate raw_comp data
  raw_data = list(comm, rec, comm, rec, comm_cwt, rec_cwt, comm_cwt, rec_cwt)
) %>% 
  mutate(
    grouping_col = case_when(
      grouping == "pst" ~ "pst_agg",
      grouping == "can" ~ "reg1"
    ),
    data = pmap(list(sample, grouping_col, raw_data), .f = clean_comp)
  )

## PREP INPUTS ----------------------------------------------------------------

prep_multinomial_inputs <- function(comp_in, 
                     # month_range = c(1, 12), 
                     data_type) {
  # gsi_trim <- gsi_in %>% 
  #   filter(!month_n < month_range[1],
  #          !month_n > month_range[2]) %>% 
  #   droplevels() %>% 
  #   dplyr::select(id, region, area, year, month, season, agg, agg_prob, pres)
  
  gsi_trim <- comp_in %>% 
    group_by(region, month) %>%
    mutate(nn = length(unique(id))) %>%
    filter(!nn < 100) %>%
    ungroup() %>%
    droplevels() %>%
    select(id, region, area, year, month, season, agg, #agg_prob, 
           pres) 
  
  # dummy dataset to replace missing values 
  dum <- gsi_trim %>% 
    group_by(region) %>% 
    nest() %>%
    droplevels() %>% 
    transmute(dum_data = map(data, function(x) {
      expand.grid(month = unique(x$month),
                  agg = unique(x$agg),
                  pres = 1)
    })
    ) %>% 
    unnest(cols = c(dum_data))
  
  # dum <- expand.grid(
  #   month = unique(gsi_trim$month),
  #   region = unique(gsi_trim$region),
  #   agg = unique(gsi_trim$agg),
  #   pres = 1)
  rand_yrs <- sample(unique(gsi_trim$year), size = nrow(dum), replace = TRUE)
  dum$year <- rand_yrs
  
  suppressWarnings(
    gsi_trim_no0 <- gsi_trim %>%
      full_join(., dum, by = c("year", "month", "agg", "pres", "region")) %>%
      mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
      droplevels()
    ) 
  
  raw_freq_table <- table(gsi_trim$agg, gsi_trim$month, gsi_trim$region)
  freq_table <- table(gsi_trim_no0$agg, gsi_trim_no0$month, gsi_trim_no0$region)
  
  gsi_wide <- gsi_trim_no0 %>% 
    pivot_wider(., names_from = agg, values_from = pres) %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
  y_obs <- gsi_wide %>% 
    select(-c(id:dummy_id)) %>% 
    as.matrix()
  head(y_obs)
  
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
  parameters <- list(z_rfac = rep(0, times = length(unique(yr_vec))),
                     z_ints = matrix(0, nrow = ncol(fix_mm), 
                                     ncol = ncol(y_obs) - 1),
                     log_sigma_rfac = 0)
  
  list("fix_mm" = fix_mm, "pred_dat" = pred_dat, 
       # "fac_key" = fac_key,
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "long_data" = gsi_trim, "long_data_no0" = gsi_trim_no0, 
       "wide_data" = gsi_wide, "freq_table" = freq_table, 
       "raw_freq_table" = raw_freq_table)
}

# prep different datasets then check catch breakdown
# NOTE: insufficient samples to include sublegal 
# pst_comm_trim <- prep_gsi(comm_pst, data_type = "comm_pst") 
# pst_legal_trim <- prep_gsi(rec_pst %>% 
#                          filter(legal == "legal"), 
#                          data_type = "legal_pst")
# can_comm_trim <- prep_gsi(comm_can, data_type = "comm_can") 
# can_legal_trim <- prep_gsi(rec_can %>% 
#                              filter(legal == "legal"), 
#                            data_type = "legal_can")
# gsi_list <- list("pst_comm" = pst_comm_trim, "pst_rec" = pst_legal_trim, 
#                  "can_comm" = can_comm_trim, "can_rec" = can_legal_trim)


# FIT --------------------------------------------------------------------------

compile(here::here("src", "multinomial_hier2.cpp"))
dyn.load(dynlib(here::here("src", "multinomial_hier2")))

fit_model <- function(x) {
  ## Make a function object
  # x <- comp2$model_inputs[[1]]
  obj <- MakeADFun(data = x$data, 
                   parameters = x$parameters, 
                   random = c("z_rfac"), 
                   DLL = "multinomial_hier2")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  f_name <- paste(x$data_type, "multinomial_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generated_data", "model_fits", f_name))
  
  return(ssdr)
}


# add model parameters and fit (exclude CWT for now)
comp2 <- comp %>%
  filter(sample == "gsi") %>% 
  mutate(model_inputs = map2(data, dataset,
                             .f = prep_multinomial_inputs),
         ssdr = map(model_inputs, .f = fit_model))

tt <- fit_model(comp2$model_inputs[[4]])

## PLOT PREDICTIONS ------------------------------------------------------------

ssdr_list <- map(comp2$model_inputs, function (x) {
  f_name <- paste(x$data_type, "multinomial_ssdr.RDS", sep = "_")
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
  pred_mat <- ssdr[rownames(ssdr) %in% "pred_probs", ]
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
    summarize(samp_g = length(unique(id))) %>% 
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
  
  f_name <- paste(x$data_type, "multinomial_pred.pdf", sep = "_")
  pdf(here::here("figs", "model_pred", "multinomial_only", f_name))
  print(pred_plot_raw)
  dev.off()
  
  f_name2 <- paste(x$data_type, "multinomial_pred.rds", sep = "_")
  saveRDS(pred_plot,
          here::here("figs", "model_pred", "multinomial_only", f_name2))
})


## HISTOGRAMS SHOWING GSI DIST -------------------------------------------------

# as above but without removing values below threshold
calc_max_prob2 <- function(grouped_data, dataset) {
  grouped_data %>% 
    summarize(agg_prob = sum(adj_prob)) %>% 
    arrange(id, desc(agg_prob)) %>%
    ungroup() %>% 
    group_by(id) %>% 
    mutate(max_assignment = max(agg_prob)) %>% 
    filter(!agg_prob < max_assignment) %>% 
    ungroup() %>% 
    distinct() %>% 
    mutate(dataset = dataset)
  }

rec_pst_h <- rec %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob2(., dataset = "rec_pst")
comm_pst_h <- comm %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob2(., dataset = "comm_pst")
#aggregate with Canadian CU focus
rec_can_h <- rec %>% 
  pool_aggs() %>% 
  group_by(id, reg1) %>%
  calc_max_prob2(., dataset = "rec_can")
comm_can_h <- comm %>% 
  pool_aggs() %>% 
  group_by(id, reg1) %>%
  calc_max_prob2(., dataset = "comm_can")

list(rec_pst_h, comm_pst_h, rec_can_h, comm_can_h) %>% 
  bind_rows() %>% 
  ggplot(.) +
  geom_histogram(aes(max_assignment)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~ dataset)
