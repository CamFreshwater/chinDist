## Dirichlet model fit
# June 2, 2020
# Fit dirichlet model to composition data from WCVI troll and southern BC
# rec fisheries

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
    ungroup() %>% 
    distinct() %>% 
    left_join(.,
              full_data %>% 
                select(id:area_n) %>% 
                distinct(),
              by = "id") %>% 
    arrange(id, desc(agg_prob))
    
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
      )
    )
}

rec_pst <- rec %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., rec)

comm_pst <- comm %>% 
  pool_aggs() %>% 
  group_by(id, pst_agg) %>%
  calc_max_prob(., comm)


# Trim based on dataset
trim_gen <- function(dat, month_range = c(1, 12)) {
  dat %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(id, region, area, year, month, season, agg, agg_prob) 
}

comm_trim <- trim_gen(comm_pst, month_range = c(4, 10))
legal_trim <- trim_gen(rec_pst %>% filter(legal == "legal"), 
                       month_range = c(6, 9))
sub_trim <- trim_gen(rec_pst %>% filter(legal == "sublegal"),
                     month_range = c(6, 9))
gsi_list <- list("comm" = comm_trim, "legal" = legal_trim, 
                 "sublegal" = sub_trim)

#spread and replace zero values
gsi_wide_list <- map(gsi_list, function(x) {
  x %>%
    pivot_wider(., names_from = agg, values_from = agg_prob) %>%
    mutate_if(is.numeric, ~replace_na(., 0)) %>%
    droplevels()  
})


## DATA PREP -------------------------------------------------------------------

# Prep data to pass to model
prep_gsi <- function(gsi_in, month_range = c(1, 12), data_type) {
  gsi_trim <- gsi_in %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(id, region, area, year, month, season, agg, agg_prob)
  
  gsi_wide <- gsi_trim %>%
    pivot_wider(., names_from = agg, values_from = agg_prob) %>%
    mutate_if(is.numeric, ~replace_na(., 0)) %>%
    droplevels() 
  
  y_obs <- gsi_wide %>% 
    select(-c(id:season)) %>% 
    as.matrix() %>% 
    Compositional::zeroreplace(., a = 1/5) 

  fix_mm <- model.matrix(~ region + month, gsi_wide) #fixed covariates only
  yr_vec <- as.numeric(gsi_wide$year) - 1
  
  # model matrix for predictions
  pred_dat <- gsi_wide %>%
    select(region, month) %>%
    distinct() %>%
    arrange(region, month)
  pred_mm <- model.matrix(~ region + month, pred_dat)
  
  data <- list(y_obs = y_obs, #obs
               rfac = yr_vec, #random intercepts
               fx_cov = fix_mm, #fixed cov model matrix
               pred_cov = pred_mm
  ) 
  
  # prep initial par values
  P <- ncol(fix_mm) #number of fixed parameters
  K <- ncol(y_obs) - 1  #number of groups - 1
  beta_in <- matrix(rnorm(n = P * K, mean = 0), P, K)
  rfac_in <- rnorm(length(unique(yr_vec)), 0 , 1)
  parameters <- list(z_ints = beta_in,
                     log_phi = runif(1, 1, 10),
                     z_rfac = rfac_in,
                     log_sigma_rfac = runif(1, 1, 10)
  )
  
  list("fix_mm" = fix_mm, "pred_mm" = pred_mm, 
       "data" = data, "parameters" = parameters, "data_type" = data_type,
       "long_data" = gsi_trim, "wide_data" = gsi_wide)
}

comm_trim <- prep_gsi(comm_pst, month_range = c(4, 10), data_type = "comm")
legal_trim <- prep_gsi(rec_pst %>% filter(legal == "legal"), 
                       month_range = c(6, 9), data_type = "legal")
sub_trim <- prep_gsi(rec_pst %>% filter(legal == "sublegal"),
                     month_range = c(6, 9), data_type = "sublegal")
gsi_list <- list("comm" = comm_trim, "legal" = legal_trim, 
                 "sublegal" = sub_trim)



## FIT MODEL -------------------------------------------------------------------

compile(here::here("src", "dirichlet_randInt_v2.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_randInt_v2")))



## Make a function object
obj <- MakeADFun(data, parameters, random = c("z_rfac"), 
                 DLL = "dirichlet_randInt_v2")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr


## PREDICTIONS -----------------------------------------------------------------

k <- ncol(y_obs) # number of stocks
stk_names <- colnames(y_obs)
N <- nrow(y_obs)

pred_mat <- ssdr[rownames(ssdr) %in% "pred_est", ]
pred_ci <- data.frame(pred_est = ssdr[rownames(ssdr) %in% "pred_est", "Estimate"],
                      pred_se =  ssdr[rownames(ssdr) %in% "pred_est", "Std. Error"],
                      region = pred_dat$region,
                      month = pred_dat$month) %>% 
  mutate(pred_prob_low = pred_est + (qnorm(0.025) * pred_se),
         pred_prob_up = pred_est + (qnorm(0.975) * pred_se)) %>%
  glimpse()

raw_prop <- gsi_trim %>% 
  group_by(region, month, year) %>%
  mutate(n_ind = length(unique(id))) %>% 
  ungroup() %>%
  group_by(region, month, year, agg, n_ind) %>%
  mutate(sum_prob = sum(agg_prob),
         mean_prob = sum_prob / n_ind) %>% 
  ungroup() %>% 
  select(region, month, year, agg, n_ind, sum_prob, mean_prob) %>% 
  distinct()

ggplot() +
  geom_point(data = raw_prop, 
             aes(x = month, y = mean_prob, fill = region),
             shape = 21, alpha = 0.4, position = position_dodge(0.6)) +
  # geom_pointrange(data = aes(x = as.factor(month), y = pred_prob, ymin = pred_prob_low,
  #                     ymax = pred_prob_up, fill = region), 
  #                 shape = 21, size = 0.4) +
  labs(y = "Probability", x = "Month") +
  facet_wrap(~ agg) +
  ggsidekick::theme_sleek()
  










# Prep data to pass to model
y_obs <- gsi_wide %>% 
  select(-c(id:season)) %>% 
  as.matrix() %>% 
  Compositional::zeroreplace(., a = 1/5) 
head(y_obs)

fix_mm <- model.matrix(~ region + month, gsi_wide) #fixed covariates only
yr_vec <- as.numeric(gsi_wide$year) - 1

# model matrix for predictions
pred_dat <- gsi_wide %>%
  select(region, month) %>%
  distinct() %>%
  arrange(region, month)
pred_mm <- model.matrix(~ region + month, pred_dat)

data <- list(y_obs = y_obs, #obs
             rfac = yr_vec, #random intercepts
             fx_cov = fix_mm, #fixed cov model matrix
             pred_cov = pred_mm
) 

# prep initial par values
P <- ncol(fix_mm) #number of fixed parameters
K <- ncol(y_obs) - 1  #number of groups - 1
beta_in <- matrix(rnorm(n = P * K, mean = 0), P, K)
rfac_in <- rnorm(length(unique(yr_vec)), 0 , 1)
parameters <- list(z_ints = beta_in,
                   log_phi = runif(1, 1, 10),
                   z_rfac = rfac_in,
                   log_sigma_rfac = runif(1, 1, 10)
)