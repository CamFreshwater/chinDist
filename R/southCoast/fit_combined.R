## Combined model fits
# April 3, 2020
# Fit combined multionomial/tweedie or multinomial/nb model to 
# stock composition and abundance data
# aggregate probability reflects summed probabilities of a given region of 
# origin (ie reg1 or 3) for a given individual

library(tidyverse)
library(TMB)
library(ggplot2)


# Import Catch -----------------------------------------------------------------
# comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
#                           "reg1RollUpCatchProb_FraserB.RDS"))
# comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
#                           "reg3RollUpCatchProb.RDS"))
comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                               "pstAggRollUpCatchProb.RDS")) %>%
  rename(regName = pstName) %>%
  filter(!regName == "NBC_SEAK")

# month range dictated by ecological scale
month_range = seq(1, 10, by = 1)

# pull months that are common to both strata and subset
# comm_months <- comp %>% 
#   select(catchReg, month) %>% 
#   distinct() %>% 
#   split(., .$catchReg) %>% 
#   map(., function(x) x %>% pull(as.numeric(month))) %>% 
#   Reduce(intersect, .)

#add dummy catch data for one month that's missing logbook data based on observed
# catch in comp data
min_catch <- comp %>% 
  filter(catchReg == "SWVI",
         month == "7") %>%
  group_by(year) %>% 
  summarise(n = length(unique(id)),
            #scalar (3) represents est effort given average CPUE
            n_days = length(unique(jDay)) * 3)
min_catch_dat <- data.frame(catch = min_catch$n,
                            catchReg = "SWVI",
                            area = NA,
                            month = "7",
                            jDay = NA,
                            year = min_catch$year,
                            boatDays = min_catch$n_days) %>% 
  mutate(cpue = catch / boatDays) %>% 
  select(catchReg:year, catch, boatDays, cpue)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                            "dailyCatch_WCVI.rds")) %>% 
  rbind(., min_catch_dat) %>% 
  mutate(reg = factor(catchReg),
         area = as.factor(area),
         month_n = as.numeric(month),
         month = as.factor(month_n),
         year = as.factor(year)
  ) %>% 
  filter(!is.na(cpue),
         #first constrain by range
         month_n %in% month_range
         # ,
         # #then drop missing months
         # month_n %in% comm_months
         ) %>% 
  droplevels() %>% 
  mutate(eff_z = as.numeric(scale(boatDays)),
         eff_z2 = eff_z^2
  ) %>% 
  arrange(catchReg, area, month)


# Clean Genetics and Prep Inputs -----------------------------------------------
source(here::here("R", "functions", "clean_composition_dat.R"))
comp_wide_l <- clean_comp(comp, month_range = month_range, check_tables = T)
comp_wide_l$tables
comp_wide <- comp_wide_l$data

fac_dat <- comp_wide %>% 
  mutate(facs = as.factor(paste(as.character(catchReg), 
                                as.character(month), sep = "_")),
         facs = fct_relevel(facs, "NWVI_1", "NWVI_2", "NWVI_3", "NWVI_4", 
                            "NWVI_5", "NWVI_6",  "NWVI_7", "NWVI_8", "NWVI_9", 
                            "NWVI_10", "NWVI_11", "NWVI_12", "SWVI_1",  
                            "SWVI_2", "SWVI_3",  "SWVI_4", "SWVI_5", "SWVI_6",  
                            "SWVI_7", "SWVI_8",  "SWVI_9", "SWVI_10", 
                            "SWVI_11", "SWVI_12"),
         facs_n = as.numeric(facs) - 1) %>% 
  select(catchReg, month, facs, facs_n)
fac_key_eff <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n) %>% 
  mutate(eff_z = mean(catch$eff_z),
         eff_z2 = mean(catch$eff_z2))

source(here::here("R", "functions", "prep_tmb_dat.R"))
inputs <- tmb_dat(catch, comp_wide, fac_dat, fac_key_eff, mod = "nb")
#check predictive matrix
head(inputs$data$X1_ij)
head(inputs$data$X1_pred_ij)


# Fit Model --------------------------------------------------------------------

## Make a function object
# compile(here::here("R", "southCoast", "tmb", "tweedie_multinomial_1re.cpp"))
# dyn.load(dynlib(here::here("R", "southCoast", "tmb", 
#                            "tweedie_multinomial_1re")))
# obj <- MakeADFun(data, parameters, random = c("z1_k", "z2_k"), 
#                  DLL = "tweedie_multinomial_1re")

compile(here::here("R", "southCoast", "tmb", "nb_multinomial_1re.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", 
                           "nb_multinomial_1re")))

obj <- MakeADFun(inputs$data, inputs$parameters, random = c("z1_k", "z2_k"), 
                 DLL = "nb_multinomial_1re")

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr

# saveRDS(ssdr, here::here("generatedData", "model_fits", "nbmult_ssdr_pst.RDS"))
# saveRDS(ssdr, here::here("generatedData", "model_fits", "nbmult_ssdr_frB.RDS"))


# PREDICTIONS ------------------------------------------------------------------

#frB versions have a different subset of stocks and months at finer spatial scale
ssdr <- readRDS(here::here("generatedData", "model_fits", "nbmult_ssdr_pst.RDS"))
# ssdr <- readRDS(here::here("generatedData", "model_fits", "twmult_ssdr_agg.RDS"))

log_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] #log pred of abundance
logit_probs <- ssdr[rownames(ssdr) %in% "logit_pred_prob", ] #logit probs of each category
pred_abund <- ssdr[rownames(ssdr) %in% "pred_abund_mg", ] #pred abundance of each category

comp_trim <- comp_wide_l$long_data
n_groups <- length(unique(comp_trim$regName))

pred_ci <- data.frame(stock = as.character(rep(unique(comp_trim$regName),
                                               each = 
                                                 length(unique(fac_key_eff$facs_n)))), 
                      logit_prob_est = logit_probs[ , "Estimate"],
                      logit_prob_se =  logit_probs[ , "Std. Error"]) %>%
  mutate(facs_n = rep(fac_key_eff$facs_n, times = n_groups)) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_low = plogis(logit_prob_est +
                                  (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est +
                                 (qnorm(0.975) * logit_prob_se)),
         abund_est = pred_abund[ , "Estimate"],
         abund_se =  pred_abund[ , "Std. Error"],
         abund_low = abund_est + (qnorm(0.025) * abund_se),
         abund_up = abund_est + (qnorm(0.975) * abund_se)) %>%
  left_join(., fac_key_eff, by = c("facs_n")) %>% 
  mutate(stock = fct_reorder(stock, desc(pred_prob)))

# calculate raw summary data for comparison
raw_prop <- comp %>% 
  filter(month_n %in% month_range) %>% 
  group_by(catchReg, month, year, regName) %>%
  summarize(samp_g = length(unique(id))) %>% 
  group_by(catchReg, month, year) %>%
  mutate(samp_total = sum(samp_g)) %>% 
  ungroup() %>% 
  mutate(samp_g_ppn = samp_g / samp_total,
         stock = fct_reorder(regName, desc(samp_g_ppn))) 

raw_abund <- catch %>% 
  group_by(catchReg, month, year) %>%
  summarize(sum_catch = sum(catch),
            sum_effort = sum(boatDays),
            agg_cpue = sum_catch / sum_effort) %>% 
  ungroup() %>% 
  mutate(month = as.character(month)) %>% 
  left_join(., raw_prop, by = c("catchReg", "month", "year")) %>% 
  mutate(catch_g = samp_g_ppn * sum_catch,
         cpue_g = samp_g_ppn * agg_cpue,
         catchReg = as.factor(catchReg), 
         month = fct_relevel(as.factor(month), "10", after = Inf)) %>% 
  filter(!is.na(stock))


# combined estimates of stock-specific CPUE
# note that effort is standardized differently between raw data and predictions,
# not apples to apples comparison
ggplot() +
  geom_point(data = raw_abund, aes(x = month, y = cpue_g, fill = catchReg),
             shape = 21, alpha = 0.3, position = position_dodge(0.6)) +
  geom_pointrange(data = pred_ci, 
                  aes(x = month, y = abund_est, ymin = abund_low, 
                      ymax = abund_up, fill = catchReg), 
                  shape = 21, position = position_dodge(0.6)) +
  facet_wrap(~stock, ncol = 2, scales = "free_y") +
  scale_fill_viridis_d(option = "C", name = "Catch Region") +
  labs(x = "Month", y = "Predicted Catch") +
  ggsidekick::theme_sleek() +
  theme(legend.position="top")

# estimates of stock compostion
ggplot() +
  geom_point(data = raw_prop, 
             aes(x = month, y = samp_g_ppn, fill = catchReg),
             shape = 21, alpha = 0.3, position = position_dodge(0.6)) +
  geom_pointrange(data = pred_ci, 
                  aes(x = month, y = pred_prob, ymin = pred_prob_low, 
                      ymax = pred_prob_up, fill = catchReg), 
                  shape = 21, position = position_dodge(0.6)) +
  facet_wrap(~stock, ncol = 2, scales = "free_y") +
  scale_fill_viridis_d(option = "C", name = "Catch Region") +
  labs(x = "Month", y = "Predicted Encounter Probability") +
  ggsidekick::theme_sleek() +
  theme(legend.position="top")

# estimates of aggregate CPUE (replace with simulated approach below)
# agg_abund <- data.frame(
#   facs_n = fac_key_eff$facs_n,
#   raw_abund_est = log_pred[ , "Estimate"],
#   raw_abund_se = log_pred[ , "Std. Error"]) %>% 
#   mutate(
#     raw_mu = exp(raw_abund_est),
#     raw_abund_low = exp(raw_abund_est + (qnorm(0.025) * raw_abund_se)),
#     raw_abund_up = exp(raw_abund_est + (qnorm(0.975) * raw_abund_se))
#   ) %>% 
#   left_join(pred_ci, ., by = "facs_n")

# catch <- catch %>% 
#   group_by(month, catchReg) %>% 
#   mutate(month_eff = mean(boatDays)) %>% 
#   ungroup() %>% 
#   mutate(catch_z = catch / month_eff)
# 
# ggplot() +
#   geom_boxplot(data = catch %>% filter(!catch == 0), aes(x = month, y = cpue),  
#              alpha = 0.4) +
#   geom_pointrange(data = agg_abund, aes(x =  month, y = raw_mu,
#                                   ymin = raw_abund_low,
#                                   ymax = raw_abund_up), color= "red") +
#   facet_wrap(~ catchReg, nrow = 2, scales = "free_y") +
#   ggsidekick::theme_sleek()


## generate better comparison of abundance model predictions by using betas
## to calculate catch across range of effort values

#betas for abundance model
abund_b <- ssdr[rownames(ssdr) %in% "b1_j", ]
pred_catch <- fac_key_eff %>% 
  mutate(
    #add model estimates
    int = abund_b[1],
    reg_b = case_when(
      catchReg == "NWVI" ~ 0,
      catchReg == "SWVI" ~ abund_b[2]
    ),
    month_b = rep(c(0, abund_b[3:(nrow(abund_b) - 2)]), times = 2),
    eff_b = abund_b[nrow(abund_b) - 1],
    eff2_b = abund_b[nrow(abund_b)]) %>% 
  split(., .$facs_n) %>% 
  map(., function(x) {
    mm <- catch %>% 
      filter(month == x$month) 
    expand_grid(x, z_eff = sample(mm$eff_z, size = 50, replace = T))
  }) %>% 
  bind_rows() %>% 
  #add estimates
  mutate(z_eff2 = z_eff^2,
         log_catch = int + reg_b + month_b + (eff_b * z_eff) + 
           (eff2_b * z_eff2),
         catch = exp(log_catch),
         dataset = "pred")

catch %>% 
  mutate(dataset = "obs") %>% 
  select(catch, dataset, catchReg, month) %>% 
  rbind(., 
        pred_catch %>% 
          select(catch, dataset, catchReg, month)
        ) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = month, y = catch, fill = dataset)) +
  facet_wrap(~ catchReg, nrow = 2, scales = "free_y") +
  ggsidekick::theme_sleek() +
  labs(fill = "Data")
# looks good! some deviations because effort isn't stratified by region and 
# July data are wonky, but otherwise solid

## export plotting data for Rmd 
list(catch = catch, pred_ci = pred_ci, raw_prop = raw_prop, raw_abund = raw_abund,
     pred_catch = pred_catch) %>%
  saveRDS(., here::here("generatedData", "model_fits", "pst_plot_list.RDS"))


## estimates of effort effects on catch
n_betas <- length(coef(m1))
eff_b <- abund_b[c(1, n_betas - 1, n_betas) , 1]

ggplot(catch) +
  geom_point(aes(x = eff_z, y = catch), 
             alpha = 0.2) +
  # lims(x = c(0, 4)) +
  stat_function(fun = function(x) exp(eff_b[1] + eff_b[2]*x + eff_b[3]*x^2)) 



# look at predicted cumulative abundance
plot_cum_dens <- function(dat, station = "JDF") {
  dat %>% 
    filter(project_name == station) %>% 
    ggplot(., aes(yday, colour = CU)) +
    stat_ecdf()
}

pred_ci %>% 
  arrange(month) %>% 
  group_by(stock, catchReg) %>% 
  mutate(total_abund = sum(abund_est),
         cum_abund = cumsum(abund_est) / total_abund) %>%
  ungroup() %>% 
  ggplot(., aes(x = month, y = cum_abund, colour = stock)) +
  geom_line() +
  geom_point() +
  facet_wrap(~catchReg)

# temp subset of stocks for comparison with cwt

stk_subset <- c("FR-early", "FR-late", "PSD", "CR-tule", "CR-bright", "WCVI")
plot_list <- map(list(pred_ci, raw_prop, raw_abund), function(x) 
  x %>% filter(stock %in% stk_subset)
)

prop_plot <- ggplot() +
  geom_pointrange(data = plot_list[[1]], aes(x = month, y = pred_prob, 
                                      ymin = pred_prob_low,
                                      ymax = pred_prob_up),
                  col = "red") +
  geom_point(data = plot_list[[2]],
             aes(x = month, y = samp_g_ppn),
             alpha = 0.4) +
  facet_wrap(stock ~ catchReg, nrow = n_groups, scales = "free_y") +
  ggsidekick::theme_sleek()

pdf(here::here("figs", "model_pred", "gsi_nb_prop_pred.pdf"))
prop_plot
dev.off()
