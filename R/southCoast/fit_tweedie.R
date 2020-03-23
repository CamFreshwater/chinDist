## Tweedie model fit
# March 20, 2020
# Fit gamma model to CPUE data from WCVI troll fishery

library(tidyverse)
library(TMB)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds")) %>% 
  filter(!is.na(cpue)) %>% 
  mutate(reg = factor(catchReg),
         area = as.factor(area),
         month = as.factor(month),
         month_n = as.numeric(month),
         year = as.factor(year),
         # z_eff = as.numeric(scale(boatDays))
         pos_catch = case_when(
           catch > 0 ~ 1,
           catch == 0 ~ 0
         )) %>% 
  rename(eff = boatDays)
catch_t <- catch %>% filter(pos_catch == 1)

fct_to_tmb_num <- function(x) {
  as.numeric(as.factor(as.character(x))) - 1
}


N <- 100
fix_dat <- data.frame(reg = as.factor(sample(c(1, 2), size = N, replace = T)),
                      month = as.factor(sample(c(1, 2), size = N, replace = T)))

# fixed intercepts
reg_ints <- c(0.25, 1) #i.e. how region 2 differs from reference 
reg_ints <- c(0.25, 1) #i.e. how region 2 differs from reference 

# model matrix for fixed effects
fix_mm <- model.matrix(~ reg : month - 1, fix_dat)



# Prep data to pass to model
y_obs <- log(catch$cpue)
yr_vec <- as.numeric(gsi_wide$year) - 1
fix_mm <- model.matrix(~ statArea + month, gsi_wide) #fixed covariates only

#make combined factor levels (necessary for increasing speed of prob. estimates)
fac_dat <- gsi_wide %>% 
  mutate(facs = as.factor(paste(statArea, month, year, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
  select(statArea, month, year, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)
