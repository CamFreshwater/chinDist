## Explore catch data from Hargreaves sampling in Barkley Sound
# August 8, 2019

library(tidyverse)
library(ggplot2)
library(ggmap)

sets <- read.csv(here::here("data", "hargreaves", "barkleySets_cleaned.csv"), 
                 stringsAsFactors = FALSE)
fish <- read.csv(here::here("data", "hargreaves", "barkleyFish_cleaned.csv"), 
                 stringsAsFactors = FALSE)
cwt <- read.csv(here::here("data", "hargreaves", "barkleyCWT_cleaned.csv"), 
                 stringsAsFactors = FALSE) %>% 
  filter(!is.na(hatchery)) %>% 
  mutate(marineRes = case_when(
    !is.na(minRelJulian) ~ julianDay - ((maxRelJulian + minRelJulian) / 2),
    maxRelJulian < julianDay ~ as.numeric(julianDay - maxRelJulian),
    maxRelJulian >= julianDay ~ NA_real_
  ))


## Plot catch locations for CK
#disaggregate by counts so that density can be plotted
mapSets <- sets %>% 
  mutate(ids = map(count, seq_len)) %>% 
  unnest() %>% 
  select(id, lat:day, ids) 

# map specs
nAm <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))
longRange = c(-126, -124.5)
latRange = c(48.75, 49.5)

blankMap <- ggplot() +
  geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
           color = "black", fill = "gray80") + 
  lims(x = longRange, y = latRange) +
  samSim::theme_sleekX()

blankMap +
  stat_density2d(data = mapSets, aes(long, lat, fill = ..level..),
                 geom = "polygon", alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) + 
  facet_wrap(~year)

blankMap +
  geom_point(data = fish, aes(long, lat, color = CWT)) + 
  facet_wrap(~year)

blankMap +
  geom_point(data = fish, aes(long, lat, color = as.factor(month))) + 
  facet_wrap(~year)

## Proportional marking and cwt rate
dum <- fish %>% 
  group_by(year, month, CWT) %>%
  tally(name = "nn")

ggplot(dum, aes(x = as.factor(month), y = nn, fill = CWT)) +
  geom_bar(stat = "identity") +
  facet_wrap(~year)


## Size distribution
ggplot(fish, aes(x = as.factor(month), y = fl, fill = CWT)) +
  geom_boxplot() +
  facet_wrap(~year)

# add quantile size distributions
fishSize <- fish %>% 
  mutate(sizeClass = findInterval(fl, sizeRange))

blankMap +
  geom_point(data = fishSize, aes(long, lat, color = sizeClass)) + 
  scale_color_viridis_c() +
  facet_wrap(~month)

## CWT 
dum <- cwt %>% 
  filter(!hatchery == "0104-ROBERTSON CREEK") %>% 
  group_by(hatchery, month) %>%
  tally(name = "nn")

ggplot(dum, aes(x = as.factor(month), y = nn, fill = hatchery)) +
  geom_bar(stat = "identity")

ggplot(cwt %>% filter(hatchery == "0104-ROBERTSON CREEK"), 
       aes(x = as.factor(month), y = maxRelJulian)) +
  geom_boxplot() +
  facet_wrap(~year)

ggplot(cwt, aes(x = as.factor(month), y = marineRes)) +
  geom_boxplot() +
  facet_wrap(~year)

