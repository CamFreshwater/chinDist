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
  
  colnames(out)[2] <- "agg"
  
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

#large scale PST groups
rec_pst <- rec %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., rec, thresh = 0.75)
comm_pst <- comm %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., comm, thresh = 0.75)
#aggregate with Canadian CU focus
rec_can <- rec %>% 
  pool_aggs() %>% 
  group_by(id, reg1) %>%
  calc_max_prob(., rec, thresh = 0.75)
comm_can <- comm %>% 
  pool_aggs() %>% 
  group_by(id, reg1) %>%
  calc_max_prob(., comm, thresh = 0.75)


## PREP MODEL INPUTS -----------------------------------------------------------

prep_gsi <- function(gsi_in, month_range = c(1, 12), data_type) {
  gsi_trim <- gsi_in %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(id, region, area, year, month, season, agg, agg_prob, pres)
  
  # dummy dataset to replace missing values 
  dum <- expand.grid(
    month = unique(gsi_trim$month),
    region = unique(gsi_trim$region),
    agg = unique(gsi_trim$agg),
    pres = 1)
  rand_yrs <- sample(unique(gsi_trim$year), size = nrow(dum), replace = TRUE)
  dum$year <- rand_yrs
  
  gsi_trim_no0 <- gsi_trim %>%
    full_join(., dum, by = c("year", "month", "agg", "pres", "region")) %>%
    mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
    droplevels() 
  
  gsi_wide <- gsi_trim_no0 %>% 
    pivot_wider(., names_from = agg, values_from = pres) %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
  raw_freq_table <- table(gsi_trim$agg, gsi_trim$month, gsi_trim$region)
  freq_table <- table(gsi_trim_no0$agg, gsi_trim_no0$month, gsi_trim_no0$region)
  
  y_obs <- gsi_wide %>% 
    select(-c(id:dummy_id)) %>% 
    as.matrix()
  head(y_obs)
  
  yr_vec <- as.numeric(gsi_wide$year) - 1
  fix_mm <- model.matrix(~ region + month, gsi_wide) #fixed covariates only
  
  #make combined factor levels (necessary for increasing speed of prob. estimates)
  fac_dat <- gsi_wide %>% 
    mutate(facs = as.factor(paste(region, as.numeric(month), sep = "_")),
           #facs = as.factor(paste(statArea, month, year, sep = "_")),
           facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
    select(region, month, facs, facs_n)
  fac_key <- fac_dat %>% 
    distinct() %>% 
    arrange(facs_n)
  
  data <- list(y_obs = y_obs, #obs
               rfac = yr_vec, #random intercepts
               fx_cov = fix_mm, #fixed cov model matrix
               n_rfac = length(unique(yr_vec)), #number of random intercepts
               all_fac = fac_dat$facs_n, # vector of factor combinations
               fac_key = fac_key$facs_n #ordered unique factor combos in fac_vec
  ) 
  parameters <- list(z_rfac = rep(0, times = length(unique(yr_vec))),
                     z_ints = matrix(0, nrow = ncol(fix_mm), 
                                     ncol = ncol(y_obs) - 1),
                     log_sigma_rfac = 0)
  
  list("fix_mm" = fix_mm, #"pred_dat" = pred_dat, 
       "fac_key" = fac_key,
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "long_data" = gsi_trim, "long_data_no0" = gsi_trim_no0, 
       "wide_data" = gsi_wide, "freq_table" = freq_table, 
       "raw_freq_table" = raw_freq_table)
}

# prep different datasets then check catch breakdown
# NOTE: insufficient samples to include sublegal 
pst_comm_trim <- prep_gsi(comm_pst %>% 
                        filter(!month == "7"), 
                      month_range = c(4, 10), data_type = "comm_pst") 
pst_legal_trim <- prep_gsi(rec_pst %>% 
                         filter(legal == "legal"), 
                       month_range = c(6, 9), data_type = "legal_pst")
can_comm_trim <- prep_gsi(comm_can %>% 
                            filter(!month == "7"), 
                          month_range = c(4, 10), data_type = "comm_can") 
can_legal_trim <- prep_gsi(rec_can %>% 
                             filter(legal == "legal"), 
                           month_range = c(6, 9), data_type = "legal_can")
gsi_list <- list("pst_comm" = pst_comm_trim, "pst_rec" = pst_legal_trim, 
                 "can_comm" = can_comm_trim, "can_rec" = can_legal_trim)

map(gsi_list, function(x) print(x$raw_freq_table))


## RUN MODEL -------------------------------------------------------------------

compile(here::here("src", "multinomial_hier.cpp"))
dyn.load(dynlib(here::here("src", "multinomial_hier")))

ssdr_list <- map(gsi_list[3:4], function (x) {
  ## Make a function object
  obj <- MakeADFun(data = x$data, 
                   parameters = x$parameters, 
                   random = c("z_rfac"), 
                   DLL = "multinomial_hier")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)

  f_name <- paste(x$data_type, "multinomial_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generatedData", "model_fits", f_name))
  
  return(ssdr)
})


## PLOT PREDICTIONS ------------------------------------------------------------

ssdr_list <- map(gsi_list, function (x) {
  f_name <- paste(x$data_type, "multinomial_ssdr.RDS", sep = "_")
  readRDS(here::here("generatedData", "model_fits", f_name))
})


plot_list <- map2(gsi_list, ssdr_list, function(x, ssdr) {
  y_obs <- x$data$y_obs
  k <- ncol(y_obs) # number of stocks
  stk_names <- colnames(y_obs)
  N <- nrow(y_obs)
  fac_key <- x$fac_key
  
  logit_probs_mat <- ssdr[rownames(ssdr) %in% "logit_probs_out", ]
  stock_n_vec <- as.character(rep(1:k, each = length(unique(fac_key$facs_n))))
  stock_vec <- as.character(rep(stk_names, each = length(unique(fac_key$facs_n))))
  pred_ci <- data.frame(stock_n = stock_n_vec, 
                        stock = stock_vec, 
                        logit_prob_est = logit_probs_mat[ , "Estimate"],
                        logit_prob_se =  logit_probs_mat[ , "Std. Error"]) %>%
    mutate(facs_n = rep(fac_key$facs_n, times = k), 
           pred_prob = plogis(logit_prob_est),
           pred_prob_se = plogis(logit_prob_se),
           pred_prob_low = plogis(logit_prob_est +
                                    (qnorm(0.025) * logit_prob_se)),
           pred_prob_up = plogis(logit_prob_est +
                                   (qnorm(0.975) * logit_prob_se))) %>%
    left_join(., fac_key, by = "facs_n") %>%
    select(-logit_prob_est, -logit_prob_se, -facs, -facs_n) 
  
  # calculate raw proportion data for comparison
  raw_prop <- x$long_data_no0 %>% 
    group_by(region, month, year, agg) %>%
    summarize(samp_g = length(unique(id))) %>% 
    group_by(region, month, year) %>%
    mutate(samp_total = sum(samp_g)) %>% 
    ungroup() %>% 
    mutate(samp_g_ppn = samp_g / samp_total,
           stock = fct_reorder(agg, desc(samp_g_ppn))) 
  
  pred_plot <- ggplot() +
    geom_point(data = raw_prop, 
               aes(x = month, y = samp_g_ppn, fill = region),
               shape = 21, alpha = 0.4, position = position_dodge(0.6)) +
    geom_pointrange(data = pred_ci,
                    aes(x = as.factor(month), y = pred_prob, ymin = pred_prob_low,
                        ymax = pred_prob_up, fill = region), 
                    shape = 21, size = 0.4, position = position_dodge(0.6)) +
    labs(y = "Probability", x = "Month") +
    facet_wrap(~ stock) +
    ggsidekick::theme_sleek()
  
  f_name <- paste(x$data_type, "multinomial_pred.pdf", sep = "_")
  pdf(here::here("figs", "model_pred", "multinomial_only", f_name))
  print(pred_plot)
  dev.off()
})

