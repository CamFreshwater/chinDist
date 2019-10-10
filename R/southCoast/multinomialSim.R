### Simulated stock composition data
## Generate stock composition data with various structural properties and then
# recover with models 
## Oct. 4 2019

require(tidyverse)
require(rethinking)
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
library(gganimate)

rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

## First simulate five year's composition with three stocks
N <- 500 #total sample number
yrs <- seq(1, 5, by = 1)
N_yr <- max(yrs) #number of years
ppnStocks <- c(0.2, 0.4, 0.6)
K <- length(ppnStocks) #number of stocks

## Generate proportion for each year
dat <- data.frame(
  id = seq(1, N, by = 1),
  year = rep(yrs, each = N / N_yr),
  stock = NA
)
trueMat <- matrix(NA, nrow = N_yr, ncol = K + 1)
trueMat[ , 1] <- yrs
for (i in seq_along(yrs)) {
  error <- runif(K, 0.0001, 0.9999)
  # generate error
  trueMat[i, 2:4] <- round(samSim::ppnAgeErr(ppnStocks, 0.1, error), digits = 2)
  simStocks <- sample(1:3, size = 100, prob = trueMat[i, 2:4], replace = T)
  dat[dat$year == yrs[i], ]$stock <- simStocks
}

### BRMS ###
#intercept only model
m_stk = brm(stock ~ 1, data = dat, family = categorical, seed = 58393)

# with random intercepts by year
m_stk2 = brm(stock ~ (1|year), data = dat, family = categorical, chains = 4)

### RETHINKING ###
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


## Practice with brms and tidybayes

m_stk2 %>% 
  spread_draws(r_year__mu2[condition, term]) %>% 
  #condition corresponds to level adn term the estimated parameter
  head(10)


# Posterior estimates 
dat %>%
  data_grid(year) %>%
  add_fitted_draws(m_stk2, dpar = TRUE) %>%
  ggplot(aes(x = year, y = .value, color = .category)) +
  stat_pointinterval(position = position_dodge(width = .4)) +
  scale_size_continuous(guide = FALSE)
