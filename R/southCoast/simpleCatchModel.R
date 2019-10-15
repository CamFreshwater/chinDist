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
  mutate(aggCPUE = aggCatch / aggEffort)

ggplot(aes(x = jDay, y = aggCatch), data = dCatch) +
  geom_point()

mCatch = brm(aggCatch ~ jDay + jDay^2 + aggEffort, family = negbinomial,
             data = dCatch, chains = 4)
mCatch2 = brm(aggCPUE ~ jDay + jDay^2, data = dCatch, chains = 4)

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