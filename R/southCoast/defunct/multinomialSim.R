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
for (i in seq_along(yrs)) {
  # generate error
  trueMat[i, 2:4] <- round(ppnErr(ppnStocks, 0.2, K), digits = 2)
  simStocks <- sample(1:3, size = 100, prob = trueMat[i, 2:4], replace = T)
  dat[dat$year == yrs[i], ]$stock <- simStocks
}

### Fit multinomial model with varying intercepts and vague priors in BRMS
# m_stk = brm(stock ~ 1, data = dat, family = categorical, chains = 4)
m_stk2 = brm(stock ~ (1|year), data = dat, family = categorical, chains = 4)
stancode(m_stk2) #check out stan model

## Extract draws
#calculate mean within each year (i.e. overall intercept plus year-specific
#effect)
m_stk2 %>% 
  spread_draws(b_mu2_Intercept, r_year__mu2[year, ]) %>%
  median_qi(yearly_mean = b_mu2_Intercept + r_year__mu2)

#plot means 
m_stk2 %>% 
  spread_draws(b_mu2_Intercept, r_year__mu2[year, ]) %>%
  median_qi(yearly_mean = b_mu2_Intercept + r_year__mu2) %>% 
  ggplot(aes(y = year, x = yearly_mean, xmin = .lower, xmax = .upper)) +
  geom_pointintervalh()

#as above but with multiple probability levels
m_stk2 %>% 
  spread_draws(b_mu2_Intercept, r_year__mu2[year, ]) %>%
  median_qi(yearly_mean = b_mu2_Intercept + r_year__mu2,
            .width = c(.95, .5)) %>% 
  ggplot(aes(y = year, x = yearly_mean, xmin = .lower, xmax = .upper)) +
  geom_pointintervalh()

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
preds = grid %>% 
  add_predicted_draws(m_stk2) %>% 
  group_by(.prediction, year) %>% 
  tally() %>% #calculate the predicted frequency from posterior
  mutate(stockFreq = n / nIter)

ggplot() + 
  stat_pointintervalh(aes(y = year, x = .value, color = .category), data = fits,
                     .width = c(.66, .95), 
                     position = position_dodge(width = .4)) +
  geom_point(aes(y = year, x = stockFreq, pch = .prediction), data = preds,
             position = position_nudge(y = -0.2)) +
  scale_color_brewer()

## CHECK STAT_DIST_SLABH func
dat %>%
  data_grid(year) %>%
  add_fitted_draws(m_stk2, dpar = TRUE) %>%
  sample_draws(30) %>%
  ggplot(aes(y = stock)) +
  stat_dist_slabh(aes(dist = "categorical", arg1 = mu2, arg2 = mu3), 
                  slab_color = "gray65", alpha = 1/10, fill = NA
  ) +
  geom_point(aes(x = response), data = ABC, shape = 21, fill = "#9ECAE1", size = 2)

# Plot predicted curves
dat %>%
  # group_by(year) %>%
  data_grid(year) %>%
  add_fitted_draws(m_stk2) %>%
  ggplot(aes(x = year, y = stock
             # , color = ordered(cyl)
             )) +
  stat_lineribbon(aes(y = .value)) +
  # geom_point(data = dat) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")


m_mpg = brm(mpg ~ hp * cyl, data = mtcars)







### Equivalent models with Rethinking package
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
  data=datWide , chains=3 , cores=3 )

k <- as.numeric(coef(m_pois))
exp(k[1])/(exp(k[1])+exp(k[2])+exp(k[3]))
exp(k[2])/(exp(k[1])+exp(k[2])+exp(k[3]))
exp(k[3])/(exp(k[1])+exp(k[2])+exp(k[3]))

## Multinomial 1
m_mult1 <- map2stan(
  alist(
    stock ~ dcategorical( softmax(0,s2,s3) ),
    s2 <- b*2, # linear model for event type 2
    s3 <- b*3, # linear model for event type 3
    b ~ dnorm(0,5)
  ) ,
  data=list(stock=dat$stock), chains=3 , cores=3 )
precis(m_mult1)
inv_logit(coef(m_mult1))
#coef * 2 and 3 gives same results as brms
