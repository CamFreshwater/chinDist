## Neg bin model fit
# April 6, 2020
# Fit neg binomial model to catch data from WCVI troll fishery

library(tidyverse)
library(TMB)

# clean function 
clean_catch <- function(dat) {
  dat %>% 
    filter(!is.na(eff),
           !eff == "0") %>% 
    mutate(reg_f = factor(region),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff)),
           eff_z2 = eff_z^2,
           eff_z3 = eff_z^3) %>% 
    arrange(reg_f, month) 
}

#commercial catch data
comm_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  rename(eff = boatDays, region = catchReg) %>% 
  clean_catch(.) %>% 
  filter(!month == "7") %>% 
  mutate(month = droplevels(month))


#recreational catch data
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                 "monthlyCatch_rec.RDS")) %>% 
  # drop months that lack data from all regions
  rename(catch = mu_catch, eff = mu_boat_trips) %>% 
  mutate(cpue = catch / eff,
         region = abbreviate(region, minlength = 4)) %>% 
  clean_catch(.) %>% 
  filter(month_n > 4, month_n < 10) %>% 
  mutate(month = droplevels(month))
  
# visualize relative catch/effort
temp <- comm_catch %>% 
  select(month, year, reg_f, eff, catch, cpue) %>% 
  mutate(dat = "comm")
plot_catch <- rec_catch %>% 
  mutate(dat = kept_legal) %>% 
  select(month, year, reg_f, eff, catch, cpue, dat) %>% 
  rbind(., temp)

# ggplot(plot_catch) +
#   geom_histogram(aes(x = cpue)) +
#   facet_wrap(~dat, scales = "free")
# ggplot(plot_catch, aes(x = eff, y = catch)) +
#   geom_point(alpha = 0.4) +
#   geom_smooth(method = "gam") +
#   facet_wrap(~dat, scales = "free") +
#   ggsidekick::theme_sleek()


# Prep data to pass to model
prep_catch <- function(catch, data_type = NULL) {
  yr_vec <- as.numeric(as.factor(as.character(catch$year))) - 1
  
  # model matrix for fixed effects
  fix_mm <- model.matrix(~ reg_f + month + eff_z + eff_z2, catch)
  
  # Factor key of unique combinations to generate predictions
  fac_key <- catch %>%
    select(reg_f, month_n) %>%
    distinct() %>%
    mutate(month = as.factor(month_n),
           facs = paste(reg_f, month_n, sep = "_"),
           facs = fct_reorder2(facs, reg_f, 
                               desc(month_n)),
           facs_n = as.numeric(as.factor(facs)) - 1
    ) %>%
    arrange(facs_n)
  
  mm_pred <- model.matrix(~ reg_f + month, fac_key) %>% 
    cbind(.,
          eff_z = rep(0, n = nrow(.)),
          eff_z2 = rep(0, n = nrow(.)))
  
  data <- list(y1_i = catch$catch,
               X1_ij = fix_mm,
               factor1k_i = yr_vec,
               nk1 = length(unique(yr_vec)),
               X1_pred_ij = mm_pred
  )
  
  # Fit simple model to initialize tmb 
  m1 <- lm(log(catch + 0.0001) ~ reg_f + month + eff_z + eff_z2, data = catch)
  
  parameters = list(
    b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
    log_phi = log(1.5),
    z1_k = rep(0, length(unique(yr_vec))),
    log_sigma_zk1 = log(0.25)
  )
  
  if (is.null(data_type)) {
    data_type <- unique(catch$kept_legal)
  }
  
  list("fix_mm" = fix_mm, "fac_key" = fac_key, "mm_pred" = mm_pred, 
       "data" = data, "parameters" = parameters, "data_type" = data_type)
}

comm_list <- prep_catch(comm_catch, data_type = "comm")

# for recreational data prep separate inputs 
rec_list1 <- rec_catch %>% 
  split(., .$kept_legal) 
nms <- names(rec_list1)
rec_list <- rec_list1 %>%
  map(., prep_catch)

fishery_list <- list(comm_list, rec_list[[1]], rec_list[[2]], rec_list[[3]])
names(fishery_list) <- c("comm", nms)

map(fishery_list, function (x) head(x$data$y1_i))

for (i in seq_along(rec_list[[1]])) {
  print(head(rec_list[[1]][[i]]))
  print(head(comm_list[[i]]))
}

# FIT --------------------------------------------------------------------------
# Compile
compile(here::here("src", "negbin_1re.cpp"))
dyn.load(dynlib(here::here("src", "negbin_1re")))

ssdr_list <- vector(length = length(fishery_list), mode = "list")
for (i in seq_along(fishery_list)) { 
  dum <- fishery_list[[i]]
  obj <- MakeADFun(dum$data, dum$parameters, random = c("z1_k"), 
                   DLL = "negbin_1re")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  
  f_name <- paste(dum$data_type, "negbin_ssdr.RDS", sep = "_")
  saveRDS(ssdr, here::here("generatedData", "model_fits", f_name))
  
  ssdr_list[[i]] <- ssdr
}


# PREDICTIONS ------------------------------------------------------------------

ssdr <- readRDS(here::here("generatedData", "model_fits", "negbin_ssdr.RDS"))

ssdr <- ssdr_list[[1]]

## Plot predictions
log_pred_fe <- ssdr[rownames(ssdr) %in% "log_prediction", ]
pred_ci <- data.frame(log_pred_est = log_pred_fe[ , "Estimate"],
                      log_pred_se =  log_pred_fe[ , "Std. Error"]) %>%
  mutate(facs_n = fac_key$facs_n) %>% 
  mutate(log_pred_low = log_pred_est + (qnorm(0.025) * log_pred_se),
         log_pred_up = log_pred_est + (qnorm(0.975) * log_pred_se),
         pred_est = exp(log_pred_est),
         pred_se = exp(log_pred_se),
         pred_low = exp(log_pred_low),
         pred_up = exp(log_pred_up)) %>%
  left_join(., fac_key, by = "facs_n") 

ggplot() +
  geom_pointrange(data = pred_ci, aes(x = as.factor(month), y = log_pred_est,
                                      ymin = log_pred_low, 
                                      ymax = log_pred_up)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~reg_f)

real_preds <- ggplot() +
  geom_pointrange(data = pred_ci, aes(x = as.factor(month), y = pred_est,
                                      ymin = pred_low, ymax = pred_up)) +
  labs(x = "month", y = "predicted catch (mean effort)") +
  ggsidekick::theme_sleek() +
  facet_wrap(~catchReg)


## plot predictions across different levels of effort
# pull coeficients
abund_b <- ssdr[rownames(ssdr) %in% "b1_j", ]

# generate data for each factor level
pred_catch <- fac_key %>%
  #add model estimates
  mutate(
    int = abund_b[1],
    reg_b = case_when(
      reg_f == "NWVI" ~ 0,
      reg_f == "SWVI" ~ abund_b[2]
    ),
    month_b = rep(c(0, abund_b[3:(nrow(abund_b) - 2)]), times = 2),
    eff_b = abund_b[nrow(abund_b) - 1],
    eff2_b = abund_b[nrow(abund_b)]) %>%
  glimpse()
  split(., .$facs_n) %>% 
  # sample effort from original dataset
  map(., function(x) {
    mm <- catch %>% 
      filter(month == x$month) 
    expand_grid(x, z_eff = sample(mm$eff_z, size = 50, replace = T))
  }) %>% 
  bind_rows() %>% 
  glimpse()
  #add estimates
  mutate(z_eff2 = z_eff^2,
         log_catch = int + reg_b + month_b + (eff_b * z_eff) + 
           (eff2_b * z_eff2),
         catch = exp(log_catch),
         dataset = "pred")
  
mm_pred %>% 
  as.data.frame() %>% 
  mutate(id = seq(1, nrow(mm_pred), by = 1)) %>% 
  split(., .$id) %>% 
  map(., function(x) {
    mm <- catch %>% 
      filter(month == x$month) 
    expand_grid(x, z_eff = sample(mm$eff_z, size = 50, replace = T))
  })

var_effort_preds <- catch %>% 
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

pdf(here::here("figs", "model_pred", "neg_bin_predictions_comm.pdf"))
real_preds
var_effort_preds
dev.off()