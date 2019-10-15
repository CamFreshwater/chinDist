### Catch model
## Fit simple model estimating catch as a quadratic function of Julian day
## Oct. 11 2019

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

dCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                              "dailyCatch_WCVI.csv"),
                   stringsAsFactors = FALSE) %>% 
  group_by(year, jDay) %>% #ignore area and temporal effects for now
  summarize(aggCatch = sum(catch),
            aggEffort = sum(boatDays)) %>% 
  filter(!aggEffort == 0) %>% 
  mutate(aggCPUE = (aggCatch + 0.0001) / aggEffort)

ggplot(aes(x = jDay, y = log(aggCatch)), data = dCatch) +
  geom_point()
ggplot(aes(x = log(aggEffort), y = log(aggCatch)), data = dCatch) +
  geom_point()
ggplot(aes(x = jDay, y = log(aggCPUE)), data = dCatch) +
  geom_point()
ggplot(aes(x = jDay, y = aggCPUE), data = dCatch) +
  geom_point()

mCatch = brm(aggCatch ~ jDay + I(jDay^2) + aggEffort, 
             family = negbinomial(link = "log"), 
             prior = c(prior(normal(0, 5), class = Intercept),
                       prior(normal(0, 2), class = b)),
             data = dCatch, chains = 4, cores = 4, iter = 4000)
mCatch2 = brm(aggCPUE ~ jDay + I(jDay^2),
              family = Gamma(link = "log"),
              data = dCatch, 
              prior = c(prior(normal(0, 5), class = Intercept),
                        prior(normal(0, 2), class = b)),
              chains = 4, cores = 4, iter = 5000,
              control = list(max_treedepth = 15))
mCatch3 = brm(aggCPUE ~ jDay + I(jDay^2),
              family = Gamma(link = "identity"),
              data = dCatch, 
              prior = c(#prior(normal(0, 5), class = Intercept),
                        prior(normal(0, 5), class = b)),
              chains = 1, cores = 4, iter = 2000, inits = "0",
              control = list(max_treedepth = 15))

plot(marginal_effects(mCatch2), points=T)

posterior_samples(mCatch2) %>% 
  as_tibble() %>% 
  mutate(iteration = 1:nrow(.))

## Generate data following estimated coefficients
int = 1.81
b1 = 0.02
b2 = -0.00003
dSeq = seq(1, 365, length.out = 50)
estShape = 0.69

out <- NULL
for(i in seq_along(dSeq)) {
  temp <- data.frame(draw = seq(1, 100, by = 1),
                     d = dSeq[i], 
                     err = rgamma(100, estShape)) %>% 
    mutate(logC = (b1 * d) + (b2 * d^2) + int + 0,
           C = exp(logC))
  out <- rbind(out, temp)
}

ggplot(aes(x = d, y = logC), data = out) +
  geom_point()

#fits and predictions from model
daySeq <- tibble(jDay = seq(from = 0, to = 365, length.out = 50))
f_quad <- fitted(mCatch2, 
                 newdata = daySeq) %>%
  as_tibble() %>%
  bind_cols(daySeq)

p_quad <- predict(mCatch2, 
          newdata = daySeq) %>%
  as_tibble() %>%
  bind_cols(daySeq)  

ggplot(data = dCatch, 
       aes(x = jDay)) +
  geom_ribbon(data = p_quad, 
              aes(ymin = Q2.5, ymax = Q97.5),
              fill = "grey83") +
  geom_smooth(data = f_quad,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              fill = "grey70", color = "black", alpha = 1, size = 1/2) +
  geom_point(aes(y = aggCPUE),
             color = "navyblue", shape = 1, size = 1.5, alpha = 1/3) +

    coord_cartesian(xlim = range(dCatch$jDay)) 

dCatch %>% 
  data_grid(jDay = seq_range(jDay, n = 100),
            aggEffort = mean(aggEffort)) %>% 
  add_fitted_draws(mCatch) %>% 
  ggplot(aes(x = jDay, y = aggCatch)) +
  stat_lineribbon(aes(y = .value)) +
  geom_point(data = dCatch) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")

mtcars %>%
  group_by(cyl) %>%
  data_grid(hp = seq_range(hp, n = 51)) %>%
  add_fitted_draws(m_mpg) %>%
  ggplot(aes(x = hp, y = mpg, color = ordered(cyl))) +
  stat_lineribbon(aes(y = .value)) +
  geom_point(data = mtcars) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")



