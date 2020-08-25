## Combine predictions 
# March 19, 2020
# Combine multinomial and neg. binomial model predictions to generate stock 
# specific estimates of abundance across months and stat. areas

library(tidyverse)

catch <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                            "nbgam_preds.RDS"))
comp <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                           "multinomial_preds.RDS")) %>% 
  rename(area = statArea) %>% 
  mutate(month_n = as.numeric(as.character(month)),
         id = paste(stock, month, area, sep = "_")) %>%
  left_join(., catch, by = c("month_n", "area")) %>% 
  select(-c(stock_n, ests, pred_prob_low, pred_prob_up, month_n, low, high))

# for now assume 0 covariance among catch and stock composition parameters
calc_se <- function(x_mu, x_se, y_mu, y_se, covar = 0) {
  cov_mat <- matrix(c(x_se, covar, covar, y_se), nrow = 2, byrow = TRUE)
  RMark::deltamethod.special("prod", mean = c(x_mu, y_mu), cov = cov_mat,
                             ses = TRUE)
}

out <- comp %>% 
  group_by(id) %>% 
  mutate(mu_k = pred_prob * mu.pred,
         se_k = calc_se(pred_prob, pred_prob_se, mu.pred, se.pred),
         mu_k_low = mu_k + (qnorm(0.025) * se_k),
         mu_k_up = mu_k + (qnorm(0.975) * se_k)) 

out %>% 
  filter(area == "123",
         month == "1")

# plot
ggplot(out) +
  geom_point(aes(x = as.factor(month), y = mu.pred, fill = area), 
             shape = 21) +
  # geom_pointrange(aes(x = as.factor(month), y = pred_prob, ymin = pred_prob_low, 
  #                     ymax = pred_prob_up, fill = statArea), shape = 21) +
  labs(y = "Probability", x = "Month") +
  facet_wrap(~ stock)  
