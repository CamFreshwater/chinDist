### Simulated stock composition data
## Generate stock composition data with various structural properties and then
# recover with models 
## Oct. 4 2019

library(rethinking)
library(tidyverse)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)

rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

## Simulate five year's composition with three stocks
N <- 500 #total sample number
yrs <- seq(1, 5, by = 1)
N_yr <- max(yrs) #number of years
ppnStocks <- c(0.2, 0.4, 0.6) #average stock proportion
K <- length(ppnStocks) #number of stocks



## Add multivariate logistic error to frequency of proportions
ppnErr <- function(ppnVec, tau, k) {
  error <- runif(k, 0.0001, 0.9999)
  dum <- log(ppnVec) + tau * qnorm(error, 0, 1) #vectorized
  x <- log(ppnVec) + tau * qnorm(error, 0, 1) - (1 / k) * sum(dum)
  exp(x) / sum(exp(x))
}

## Generate proportions for each year
dat <- data.frame(
  id = seq(1, N, by = 1),
  year = rep(yrs, each = N / N_yr),
  stock = NA
)
trueMat <- matrix(NA, nrow = N_yr, ncol = K + 1)
trueMat[ , 1] <- yrs
for (i in seq_along(yrs)) {   # generate error
  trueMat[i, 2:4] <- round(ppnErr(ppnStocks, 0.2, K), digits = 2)
  simStocks <- sample(1:3, size = 100, prob = trueMat[i, 2:4], replace = T)
  dat[dat$year == yrs[i], ]$stock <- simStocks
}

### Fit multinomial model with varying intercepts and vague priors in BRMS
m_stk = brm(stock ~ 1, data = dat, family = categorical, chains = 4)
m_stk2 = brm(stock ~ (1|year), data = dat, family = categorical, chains = 4)
stancode(m_stk2) #check out stan model

#Generate posterior estimates
dat %>%
  data_grid(year) %>%
  add_fitted_draws(m_stk2) %>%
  ggplot(aes(x = year, y = .value, color = .category)) +
  stat_pointinterval(position = position_dodge(width = .4)) +
  scale_size_continuous(guide = FALSE)

#Generate posterior predictive estimates
grid = dat %>% 
  data_grid(year)
fits = grid %>% 
  add_fitted_draws(m_stk2)
nIter <- length(unique(fits$.draw))

obs <- trueMat %>% 
  as.data.frame() %>% 
  rename(year = V1) %>%
  gather(key = stock, value = freq, -year) %>% 
  mutate(stock = fct_recode(stock, "1" = "V2", "2" = "V3", "3" = "V4"))
# 
# preds = grid %>% 
#   add_predicted_draws(m_stk2) %>% 
#   group_by(.prediction, year) %>% 
#   tally() %>% #calculate the predicted frequency from posterior
#   mutate(stockFreq = n / nIter)

ggplot() + 
  stat_pointintervalh(aes(y = year, x = .value, color = .category), data = fits,
                     .width = c(.66, .95)) +
  geom_point(aes(y = year, x = freq, color = stock), data = obs, pch = 6,
             position = position_nudge(y = -0.2)) 
# +
#   scale_color_brewer()


### Equivalent models (no year effects) with Rethinking package
## Poisson multinomial 
datWide <- dat %>% 
  group_by(year, stock) %>% 
  tally() %>%
  spread(stock, n) %>% 
  dplyr::rename(stock1 = 2, stock2 = 3, stock3 = 4)

m_pois <- map2stan(
  alist(
    stock1 ~ dpois(lambda1),
    stock2 ~ dpois(lambda2),
    stock3 ~ dpois(lambda3),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    log(lambda3) <- a3,
    c(a1,a2, a3) ~ dnorm(0,100)
  ),
  data=datWide , chains=3, cores=3)
rethinking::stancode(m_pois) #check out stancode

k <- as.numeric(coef(m_pois))
exp(k[1])/(exp(k[1])+exp(k[2])+exp(k[3]))
exp(k[2])/(exp(k[1])+exp(k[2])+exp(k[3]))
exp(k[3])/(exp(k[1])+exp(k[2])+exp(k[3]))

## Multinomial 1
m_mult <- map2stan(
  alist(
    stock ~ dcategorical(softmax(0,s2,s3)),
    s2 <- b*2, # linear model for event type 2
    s3 <- b*3, # linear model for event type 3
    b ~ dnorm(0,5)
  ),
  data = list(stock = dat$stock), chains=1, cores=1)
rethinking::stancode(m_mult) #check out stancode
precis(m_mult)
#coef * 2 and 3 gives same results as brms, but parameters are difficult to 
#interpret
