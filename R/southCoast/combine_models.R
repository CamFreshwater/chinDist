## Combine predictions 
# March 19, 2020
# Combine multinomial and neg. binomial model predictions to generate stock 
# specific estimates of abundance across months and stat. areas

library(tidyverse)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                            "nbgam_preds.RDS"))
comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                           "multinomial_preds.RDS")) %>% 
  mutate(month_n = as.numeric(as.character(month))) %>%
  rename(area = statArea) %>% 
  left_join(., catch, by = c("month_n", "area")) %>% 
  select(-c(stock_n, ests, pred_prob_low, pred_prob_up, month_n, low, high))

comp %>% 
  select(pred_prob, mu.pred) %>% 
  cov()

calcSEs <- function(x, y) {
  rs_riv <- x[x$stage == "river", ]$mean
  
  x$final_mu <- NA
  x$final_se <- NA
  for (i in 1:nrow(x)) {
    if (x$stage[i] == "river") {
      x$final_mu[i] <- rs_riv
      x$final_se[i] <- x$se[i] 
    }
    if (x$stage[i] == "est" | x$stage[i] == "microtroll") {
      stg <- as.character(x$stage[i])
      cov_m <- y[c(stg, "river"), c(stg, "river")]
      
      x$final_mu[i] <- rs_riv / x[[i, "mean"]]
      x$final_se[i] <- RMark::deltamethod.special("prod",
                                                  #for vector of means pass
                                                  #untransformed value 
                                                  mean = c(1 / x[[i, "mean"]],
                                                           rs_riv),
                                                  cov = cov_m, 
                                                  ses = TRUE) 
    }
  }
  return(x)
}