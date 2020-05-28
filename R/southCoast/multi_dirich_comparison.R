## Compare dirichlet and multinomial predictions
# May 19, 2020
# Preliminary comparison using region 3 roll ups and fixed predictions only; 
# also compares TMB version to DirichletReg version

library(tidyverse)
library(TMB)
library(DirichletReg)

dat_i <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                            "reg3RollUpCatchProb.RDS")) %>% 
  filter(!month_n < 4,
         !month_n > 9) %>% 
  dplyr::select(id, statArea, year, month, month_n, regName, pres, catchReg) 
#subset probabilities to same inviduals
dat_p <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                    "reg3RollUpCatchProb_unfiltered.RDS")) %>% 
  mutate(
    month_n = as.numeric(month),
    month = as.factor(month_n),
    year =  as.factor(year),
    area_n = as.numeric(as.character(statArea)),
    catchReg = case_when(
      area_n < 125 & area_n > 27 ~ "SWVI",
      area_n < 25 ~ "SWVI",
      TRUE ~ "NWVI"
    ),
    pres = 1,
    catchReg = as.factor(catchReg), 
    statArea = as.factor(statArea)
    ) %>% 
  filter(id %in% dat_i$id,
         !month_n < 4,
         !month_n > 9) %>% 
  dplyr::select(id, statArea, year, month, month_n, regName, prob = aggProb, 
                catchReg) 


# trim to subset of months and widen
dat_p2 <- dat_p %>% 
  arrange(regName) %>% 
  pivot_wider(., names_from = regName, values_from = prob) %>% 
  mutate_if(is.numeric, ~replace_na(., 1e-5)) %>%
  droplevels()  %>% 
  arrange(id)
dat_i2 <- dat_i %>% 
  arrange(regName) %>% 
  pivot_wider(., names_from = regName, values_from = pres) %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  droplevels() %>% 
  arrange(id)


## PREP MODEL INPUTS -----------------------------------------------------------

#dirichlet observation matrix and substitute zero values
y_dir <- dat_p2 %>% 
  select(-c(id:catchReg)) %>% 
  as.matrix()
y_dir <- y_dir / rowSums(y_dir)
y_dir_scaled <- y_dir * 100

# multinomial observation matrix
y_mult <- dat_i2 %>% 
  select(-c(id:catchReg)) %>% 
  as.matrix()

# shared fixed covariate matrix
fix_mm <- model.matrix(~ (month + catchReg), dat_i2)
P <- ncol(fix_mm) - 1
K <- ncol(y_dir)


#Dirichlet TMB inputs
#initial parameter values
beta_in_d <- matrix(rnorm(n = (P + 1) * K, mean = 0), (P+1), K)

#prediction model matrix
pred_dat_d <- dat_p2 %>%
  select(month, catchReg) %>%
  distinct() %>%
  arrange(month, catchReg)
pred_mm_d <- model.matrix(~ month + catchReg, pred_dat_d)

compile(here::here("src", "dirichlet_fixInt.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_fixInt")))


# Multinomial TMB inputs
#initial parameter values
beta_in_m = matrix(rnorm(n = (P + 1) * (K - 1), mean = 0), nrow = P + 1,
                   ncol = K - 1)

fac_dat <- dat_i2 %>% 
  mutate(facs = as.factor(paste(month, catchReg, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
  select(month, catchReg, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)

compile("src/multinomial_fixInt.cpp")
dyn.load(dynlib("src/multinomial_fixInt"))


## FIT MODELS ------------------------------------------------------------------

## Fit dirichlet from TMB
obj <- MakeADFun(data=list(fx_cov = fix_mm,
                           y_obs = y_dir_scaled,
                           pred_cov = pred_mm_d),
                 parameters=list(z_ints = beta_in_d),
                 DLL="dirichlet_fixInt")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr_dir <- sdreport(obj)
ssdr_dir <- summary(sdr_dir)
betas_dr1 <- ssdr_dir[rownames(ssdr_dir) %in% "z_ints", "Estimate"]
preds_dr1 <- ssdr_dir[rownames(ssdr_dir) %in% "pred_pi_prop", "Estimate"]

#tmb predictions
pred1_mat <- matrix(as.vector(preds_dr1), nrow = nrow(pred_dat_d), ncol = K)
colnames(pred1_mat) <- colnames(y_dir_scaled)
pred_dat1 <- cbind(pred_dat_d, 
                   pred1_mat)  %>% 
  mutate(model = "tmb")


## Fit dirichlet regression from package
temp <- dat_p2 %>% 
  select(id:catchReg)
temp$comp <- DR_data(y_dir)

# fit_dr <- DirichletReg::DirichReg(comp ~ month + catchReg, data = temp)
fit_drB <- DirichletReg::DirichReg(comp ~ month + catchReg, data = temp,
                                  model = "alternative")
betas_dr2 <- coef(fit_drB) %>%
  do.call(c, .) %>%
  as.vector()

#DR predictions
pred_dat2 <- dat_p2 %>% 
  select(id:catchReg) %>% 
  cbind(., predict(fit_dr)) %>% 
  group_by(month, catchReg) %>% 
  summarize(Colmb = mean(Colmb),
            FrsrR = mean(FrsrR), 
            Other = mean(Other),
            PgtSn = mean(PgtSn)) %>% 
  ungroup() %>% 
  mutate(model = "dr")

## Fit dirichlet regression from second package
resp <- temp[, c("month", "catchReg")]
fit_dr2 <- diri.reg(y_dir, resp)
fit_dr2$be

# combine and compare
rbind(pred_dat1, pred_dat2) %>% 
  pivot_longer(., cols = Colmb:PgtSn, names_to = "stock", values_to = "prob") %>% 
  ggplot(.) +
  geom_point(aes(x = stock, y = prob, fill = model), shape = 21) +
  facet_grid(month~catchReg)


## Fit multinomial TMB with integer data 
obj <- MakeADFun(data = list(fx_cov = fix_mm,
                             y_obs = y_mult,
                             all_fac = fac_dat$facs_n,
                             fac_key = fac_key$facs_n),
                 parameters = list(z_ints = beta_in_m),
                 DLL="multinomial_fixInt")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr_m <- sdreport(obj)
ssdr_m <- summary(sdr_m)

preds_m <- ssdr_m[rownames(ssdr_m) %in% "logit_probs_out", "Estimate"]
pred2_mat <- matrix(plogis(as.vector(preds_m)), 
                    nrow = nrow(pred_dat_d), ncol = K)


## COMPARE TMB MODELS ----------------------------------------------------------
make_pred_df <- function(model, ssdr_in) {
  dum <- expand.grid(month = pred_dat_d$month,
                     catchReg = pred_dat_d$catchReg,
                     stock = levels(dat_p$regName)) %>% 
    distinct() %>%
    arrange(stock, month, catchReg)
  
  if (model == "mult") {
    out <- dum %>% 
      mutate(
        logit_prob_est = ssdr_in[ , "Estimate"],
        logit_prob_se = ssdr_in[ , "Std. Error"],
        mu = plogis(logit_prob_est),
        lo = plogis(logit_prob_est + (qnorm(0.025) * logit_prob_se)),
        up = plogis(logit_prob_est + (qnorm(0.975) * logit_prob_se)),
        model = "mult"
        # ,
        # covered = lo < true_prob & up > true_prob
      ) %>% 
      select(-logit_prob_est, -logit_prob_se)
  } 
  
  if (model == "dir") {
    out <- dum %>% 
      mutate(mu = ssdr_in[ , "Estimate"],
             se = ssdr_in[ , "Std. Error"],
             lo = mu + (qnorm(0.025) * se),
             up = mu + (qnorm(0.975) * se),
             model = "dir"
             # ,
             # covered = lo < true_prob & up > true_prob
             ) %>% 
      select(-se)
  }
  return(out)
}

dir_preds <- make_pred_df(model = "dir", 
             ssdr_in = ssdr_dir[rownames(ssdr_dir) %in% "pred_pi_prop", ])
mult_preds <- make_pred_df(model = "mult", 
             ssdr_in = ssdr_m[rownames(ssdr_m) %in% "logit_probs_out", ])

rbind(dir_preds, mult_preds) %>% 
  filter(catchReg == "NWVI") %>% 
  ggplot(.) +
  geom_pointrange(aes(x = stock, y = mu, ymin = lo, ymax = up, fill = model), 
                  shape = 21) +
  facet_wrap(~month, scales = "free_y")
