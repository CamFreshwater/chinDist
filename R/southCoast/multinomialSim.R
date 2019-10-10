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
m_stk = brm(stock ~ 1, data = dat, family = categorical, seed = 58393)
m_stk2 = brm(stock ~ (1|year), data = dat, family = categorical)



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

## brms example
dd <- data.frame(
  y1 = rbinom(N, 10, 0.3), y2 = rbinom(N, 10, 0.5), 
  y3 = rbinom(N, 10, 0.7), x = rnorm(N)
)
dd$size <- with(dd, y1 + y2 + y3)
dd$y <- with(dd, cbind(y1, y2, y3))



# # simulate career choices among 500 individuals 11.35
# N <- 500 # number of individuals
# income <- 1:3 # expected income of each career
# score <- 0.5*income # scores for each career, based on income
# # next line converts scores to probabilities
# p <- softmax(score[1],score[2],score[3])
# # now simulate choice
# # outcome career holds event type values, not counts
# career <- rep(NA,N) # empty vector of choices for each individual
# # sample chosen career for each individual
# for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )

N <- 100
# simulate family incomes for each individual
family_income <- runif(N)
# assign a unique coefficient for each type of event
b <- (1:-1)
career <- rep(NA,N) # empty vector of choices for each individual
for ( i in 1:N ) {
  score <- 0.5*(1:3) + b*family_income[i]
  p <- softmax(score[1],score[2],score[3])
  career[i] <- sample( 1:3 , size=1 , prob=p )
}

m10.17 <- map(
  alist(
    career ~ dcategorical( softmax(0,s2,s3) ),
    s2 <- a2 + b2*family_income,
    s3 <- a3 + b3*family_income,
    c(a2,a3,b2,b3) ~ dnorm(0,5)
  ) ,
  data=list(career=career,family_income=family_income) )
precis(m10.17)

inv_logit(coef(m10.17)[1:2])

## Generate implied predictions
# call link with sim data 
incSeq <- seq(from = 0, to = 1, length.out = N)

mu <- link(m10.17, 
           data = data.frame(career = career, family_income = family_income))
# summarize samples across cases
mu_mean_s2 <- apply( mu$s2 , 2 , mean )
mu_PI_s2 <- apply( mu$s2 , 2 , PI )
mu_mean_s3 <- apply( mu$s3 , 2 , mean )
mu_PI_s3 <- apply( mu$s3 , 2 , PI )

plot(career ~ family_income, col = rangi2)
lines(incSeq, mu_mean_s2)
lines(incSeq, mu_mean_s3)






library(rethinking)
data(UCBadmit)
d <- UCBadmit
# binomial model of overall admission probability
m_binom <- map(
  alist(
    admit ~ dbinom(applications,p),
    logit(p) <- a,
    a ~ dnorm(0,100)
  ),
  data=d )
# Poisson model of overall admission rate and rejection rate
d$rej <- d$reject # 'reject' is a reserved word
m_pois <- map2stan(
  alist(
    admit ~ dpois(lambda1),
    rej ~ dpois(lambda2),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    c(a1,a2) ~ dnorm(0,100)
  ),
  data=d , chains=3 , cores=3 )

logistic(coef(m_binom))
k <- as.numeric(coef(m_pois))
exp(k[1])/(exp(k[1])+exp(k[2]))





