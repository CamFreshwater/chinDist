## Combined model fits
# April 3, 2020
# Fit combined multionomial/tweedie model to stock composition and abundnace data
# aggregate probability reflects summed probabilities of a given region of 
# origin (ie reg1 or 3) for a given individual

library(tidyverse)
library(TMB)
library(ggplot2)


# Import Catch -----------------------------------------------------------------
month_range = c(4, 10) #month range (generally chosen based on gsi data)

#add dummy catch data for one month that's missing based on gsi data
min_catch <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                                "reg1RollUpCatchProb_Fraser.RDS"))  %>% 
  filter(catchReg == "SWVI",
         month == "7") %>%
  group_by(year) %>% 
  summarise(n = length(unique(flatFileID)),
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

catch <- catch_full %>% 
  rbind(., min_catch_dat) %>% 
  mutate(reg = factor(catchReg),
         area = as.factor(area),
         month_n = as.numeric(month),
         month = as.factor(month_n),
         year = as.factor(year)) %>% 
  filter(!is.na(cpue),
         !month_n < month_range[1],
         !month_n > month_range[2]) %>% 
  droplevels() %>% 
  arrange(catchReg, month)


# Import Genetics --------------------------------------------------------------
# region3 roll up
# reg3 <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
#                            "reg3RollUpCatchProb.RDS")) 
# region1 roll up
reg1_fr <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                              "reg1RollUpCatchProb_Fraser.RDS")) 


# Trim based on dataset
trim_gen <- function(dat, month_range = c(1, 12)) {
  dat %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(flatFileID, statArea, year, month, month_n, season, 
                  regName, pres, catchReg) 
}

gsi_trim <- trim_gen(reg1_fr, month_range = month_range)

table(gsi_trim$regName, gsi_trim$month)

# dummy dataset to replace missing values 
dum <- expand.grid(
  month = unique(gsi_trim$month),
  catchReg = unique(gsi_trim$catchReg),
  regName = unique(gsi_trim$regName),
  pres = 1)
#add random subset of yrs to avoid overparameterizing model
rand_yrs <- sample(unique(gsi_trim$year), size = nrow(dum), replace = TRUE)
dum$year <- rand_yrs

temp <- gsi_trim %>%
  full_join(., dum, by = c("year", "month", "regName", "pres", "catchReg")) %>%
  mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1))

# check for no zeros
table(temp$regName, temp$month, temp$catchReg)
table(catch$month, catch$catchReg)

gsi_wide <- temp %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  droplevels()


# Prep Data for TMB ------------------------------------------------------------

# fixed effects model matrices
fix_mm_c <- model.matrix(~ catchReg + month, data = catch)
fix_mm_gsi <- model.matrix(~ catchReg + month, data = gsi_wide)

# observed stock composition
y_obs <- gsi_wide %>% 
  select(-c(flatFileID:dummy_id)) %>% 
  as.matrix()
head(y_obs)

#helper function to convert factors 
fct_to_tmb_num <- function(x) {
  as.numeric(as.factor(as.character(x))) - 1
}

# fit dummy model to speed up tweedie estimates
m1 <- lm(log(catch + 0.0001) ~ catchReg + month, data = catch)

# make fixed effects factor key based on stock composition data
fac_dat <- catch %>% 
  mutate(facs = as.factor(paste(as.character(catchReg), 
                                as.character(month), sep = "_")),
         #so far unable to automate this...
         facs = fct_relevel(facs, "NWVI_10", after = 7),
         facs = fct_relevel(facs, "SWVI_10", after = Inf),
         facs_n = as.numeric(facs) - 1) %>% 
  select(catchReg, month, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)
# predictive model matrix
mm_pred <- model.matrix(~ catchReg + month, data = fac_key)

# other parameters
n_groups <- ncol(y_obs)
b1_n <- length(coef(m1))
b2_n <- ncol(mm_pred) * (n_groups - 1) #each par est for each non-ref level
# vectors of random effects
fac1k <- fct_to_tmb_num(catch$year)
fac2k <- fct_to_tmb_num(gsi_wide$year)
nk1 <- length(unique(fac1k))
nk2 <- length(unique(fac2k))

# combine
data <- list(
  #abundance data
  y1_i = catch$cpue,
  X1_ij = fix_mm_c,
  factor1k_i = fac1k,
  nk1 = nk1,
  X1_pred_ij = mm_pred,
  #composition data
  y2_ig = y_obs,
  X2_ij = fix_mm_gsi,
  factor2k_i = fac2k,
  nk2 = nk2,
  m2_all_fac = fac_dat$facs_n,
  m2_fac_key = fac_key$facs_n
)

parameters = list(
  b1_j = coef(m1) + rnorm(length(coef(m1)), 0, 0.01),
  log_phi = log(1.1),
  logit_p = boot::logit(0.8),
  z1_k = rep(0, length(unique(fac1k))),
  log_sigma_zk1 = log(0.25),
  b2_jg = matrix(0, nrow = ncol(mm_pred), ncol = (n_groups - 1)), 
  z2_k = rep(0, times = length(unique(fac2k))),
  log_sigma_zk2 = log(0.25)
)


# Fit Model --------------------------------------------------------------------

## Make a function object
compile(here::here("R", "southCoast", "tmb", "tweedie_multinomial_1re.cpp"))
dyn.load(dynlib(here::here("R", "southCoast", "tmb", 
                           "tweedie_multinomial_1re")))
obj <- MakeADFun(data, parameters, random = c("z1_k", "z2_k"), 
                 DLL = "tweedie_multinomial_1re")

opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr
ssdr <- summary(sdr)
ssdr
