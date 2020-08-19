## Dirichlet model fit
# July 17, 2020
# Fit dirichlet model to stock composition data at various aggregate levels
# Accounts for uncertain GSI and uses julian day within a strata as a sampling 
# event
# Same as fit_dirichlet but adds splines 

library(tidyverse)
library(TMB)
library(ggplot2)
library(mgcv)

# recreational data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(legal == "legal") %>%
  mutate(
    temp_strata = paste(month, region, sep = "_"),
    sample_id = paste(temp_strata, jDay, year, sep = "_"),
    min_m = case_when(
      region == "N. Strait of Georgia" ~ 5,
      region %in% c("Queen Charlotte Sound", 
                    "Queen Charlotte and\nJohnstone Straits") ~ 6,
      region == "Juan de Fuca Strait" ~ 3,
      region == "S. Strait of Georgia" ~ 1
    ),
    max_m = case_when(
      region %in% c("Juan de Fuca Strait", "S. Strait of Georgia") ~ 10,
      region == "Queen Charlotte Sound" ~ 8, 
      region %in% c("N. Strait of Georgia", 
                    "Queen Charlotte and\nJohnstone Straits")  ~ 9
    ),
    region = fct_relevel(as.factor(region), 
                         "Queen Charlotte and\nJohnstone Straits",
                         "N. Strait of Georgia", "S. Strait of Georgia", 
                         "Juan de Fuca Strait")
  ) %>% 
  group_by(region) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>% 
  filter(!region %in% c("Queen Charlotte Sound")) %>%
  ungroup() %>% 
  droplevels()

table(rec$region, rec$month)
table(rec$area, rec$month)
table(rec$year, rec$month, rec$region)

# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds")) %>%
  mutate(sample_id = paste(temp_strata, jDay, year, sep = "_"))
  

# helper function to calculate aggregate probs
calc_agg_prob <- function(grouped_data, full_data) {
  grouped_data %>% 
    summarize(agg_prob = sum(adj_prob), .groups = "drop") %>% 
    arrange(sample_id, desc(agg_prob)) %>%
    distinct() %>% 
    left_join(.,
              full_data %>% 
                #important to subset appropriately to remove individual traits
                #that lead to duplicates
                select(sample_id, region, year, month, gear, month_n) %>% 
                distinct(),
              by = "sample_id")
}

# helper function to pool non-focal stocks
pool_aggs <- function(full_data) {
  full_data %>% 
    mutate(
      # pst_agg = case_when(
      #   grepl("CR-", pst_agg) ~ "CR",
      #   grepl("CST", pst_agg) ~ "CA/OR/WA",
      #   pst_agg %in% c("NBC_SEAK", "WCVI") ~ "BC-coast",
      #   pst_agg %in% c("PSD", "SOG") ~ "SalSea",
      #   TRUE ~ pst_agg
      # ),
      reg1 = case_when(
        Region1Name %in% c("Fraser_Fall", "ECVI", "Fraser_Summer_4.1", 
                           "Fraser_Spring_4.2", "Fraser_Spring_5.2", 
                           "Fraser_Summer_5.2", "WCVI", "SOMN") ~ Region1Name,
        TRUE ~ "Other"
      )
    )
}

# helper function to use above for gsi and cwt samples (cwt ignored for now)
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

# comp_in <- comp$data[[2]]
# data_type <- comp$dataset[[2]]

prep_dir_inputs <- function(comp_in, data_type) {
  gsi_trim <- comp_in %>% 
    group_by(region, month, year) %>%
    mutate(nn = sum(agg_prob)) %>%
    #remove strata with less than 10 individuals total 
    filter(!nn < 10) %>%
    ungroup() %>%
    droplevels() %>%
    select(sample_id, region, year, month, month_n, agg, agg_prob) %>%
    distinct()
  
  gsi_wide <- gsi_trim %>% 
    pivot_wider(., names_from = agg, values_from = agg_prob) %>%
    mutate_if(is.numeric, ~replace_na(., 0.000001))
  
  y_obs <- gsi_wide %>% 
    select(-c(sample_id:month_n)) %>% 
    as.matrix() 

  yr_vec <- as.numeric(gsi_wide$year) - 1
  
  #generate model matrix based on GAM
  months <- unique(gsi_wide$month_n)
  n_months <- length(months)
  n_knots <- ifelse(max(months) == 12, 4, 3)
  spline_type <- ifelse(max(months) == 12, "cc", "tp")
  # response variable doesn't matter, since not fit
  m1 <- gam(rep(0, length.out = nrow(gsi_wide)) ~ 
              region + s(month_n, bs = spline_type, k = n_knots, by = region),
            # knots = list(month_n = c(min(months), max(months))),
            data = gsi_wide)
  fix_mm <- predict(m1, type = "lpmatrix")
  
  # data frame for predictions
  # account for strong differences in sampling months for sport fishery
  if (comp_in$gear[1] == "sport") {
    pred_dat <- split(gsi_wide, gsi_wide$region) %>% 
      map(., function (x) {
        expand.grid(
          month_n = seq(min(x$month_n), 
                        max(x$month_n),
                        by = 0.1),
          region = unique(x$region)
        )
      }) %>% 
      bind_rows()
  } else {
    pred_dat <- expand.grid(
      month_n = seq(min(gsi_wide$month_n), 
                    max(gsi_wide$month_n),
                    by = 0.1),
      region = unique(gsi_wide$region)
    )
  }
  pred_mm <- predict(m1, pred_dat, type = "lpmatrix")
  
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
  x <- comp2$model_inputs[[2]]
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
  
  f_name <- paste(x$data_type, "dirichlet_spline_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generated_data", "model_fits", f_name))
}

# fit
map(comp2$model_inputs, .f = fit_model)


## PLOT PREDICTIONS ------------------------------------------------------------

comp2$ssdr <- map(comp2$model_inputs, function (x) {
  f_name <- paste(x$data_type, "dirichlet_spline_ssdr.RDS", sep = "_")
  readRDS(here::here("generated_data", "model_fits", f_name))
})

pal <- readRDS(here::here("generated_data", "color_pal.RDS"))

plot_list <- map2(comp2$model_inputs, comp2$ssdr, function(x, ssdr) {
  x <- comp2$model_inputs[[2]]
  # ssdr <- comp2$ssdr[[2]]
  y_obs <- x$data$y_obs
  k <- ncol(y_obs) # number of stocks
  stk_names <- colnames(y_obs)
  N <- nrow(y_obs)
  pred_dat <- x$pred_dat
  
  pred_dat2 <- purrr::map_dfr(seq_len(k), ~pred_dat)
  pred_mat <- ssdr[rownames(ssdr) %in% "inv_logit_pred_pi_prop", ]
  pred_ci <- data.frame(
    stock = as.character(rep(stk_names, each = nrow(pred_dat))),
    link_prob_est = pred_mat[ , "Estimate"],
    link_prob_se =  pred_mat[ , "Std. Error"]
  ) %>% 
    cbind(pred_dat2, .) %>%
    mutate(pred_prob_est = car::logit(link_prob_est),
           pred_prob_low = pmax(0,
                                car::logit(link_prob_est + (qnorm(0.025) * 
                                                              link_prob_se))),
           pred_prob_up = car::logit(link_prob_est + (qnorm(0.975) * link_prob_se)) 
           #region = abbreviate(region, minlength = 4)
           ) 

  
  # calculate raw proportion data for comparison (wide necessary to 
  # account for infilling)
  n_stks <- length(stk_names)
  raw_prop <- x$wide_data %>% 
    pivot_longer((ncol(.) - n_stks + 1):ncol(.), 
                 names_to = "agg", values_to = "agg_prob") %>% 
    group_by(region, month, year, agg) %>%
    summarize(samp_g = sum(agg_prob)) %>% 
    group_by(region, month, year) %>%
    mutate(samp_total = sum(samp_g)) %>% 
    ungroup() %>% 
    mutate(samp_g_ppn = samp_g / samp_total,
           stock = fct_reorder(agg, desc(samp_g_ppn))#, 
           # region = abbreviate(region, minlength = 4)
           ) 
  
  # raw_prop %>%
  #   filter(region == "JndF", stock == "SOG", month == "5") %>%
  #   print(n = Inf)
  # 
  pred_plot <- ggplot(data = pred_ci, aes(x = month_n)) +
    geom_line(aes(y = pred_prob_est, colour = region )) +
    geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = region), 
                alpha = 0.5) +
    # scale_fill_manual(name = "Region", values = pal) +
    # scale_colour_manual(name = "Region", values = pal) +
    labs(y = "Predicted Encounter Probability", x = "Month") +
    scale_x_continuous(breaks = seq(1, 12, by = 1)) +
    facet_wrap(~stock, scales = "free_y") +
    ggsidekick::theme_sleek()
  
  # same as above but scales not fixed
  pred_plot2 <- pred_plot +
    facet_wrap(~stock)
  
  # same as above but with annual observations shown
  pred_plot_raw <- pred_plot +
    geom_point(data = raw_prop, 
               aes(x = as.numeric(as.character(month)), 
                   y = samp_g_ppn, fill = region),
               shape = 21, alpha = 0.4, position = position_dodge(0.6)) +
    facet_wrap(~stock)
  
  f_name <- paste(x$data_type, "dirichlet_spline_pred.pdf", sep = "_")
  pdf(here::here("figs", "model_pred", "dirichlet_only", f_name))
  print(pred_plot)
  print(pred_plot2)
  print(pred_plot_raw)
  dev.off()
  
  f_name2 <- paste(x$data_type, "dirichlet_spline_pred.rds", sep = "_")
  saveRDS(pred_plot,
          here::here("figs", "model_pred", "dirichlet_only", f_name2))
})
