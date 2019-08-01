## Initial exploration of tagging data
# July 22, 2019

library(tidyverse)
library(ggplot2)

setDat <- read.csv(here::here("data", "taggingData", "cleanSetData.csv")) 
chinDat <- read.csv(here::here("data", "taggingData", "tagChinData.csv"))

# Merge catch and fish data
calcWt <- function(length, girth) {
  #0.3937 accounts for conversion to inches; output is pounds
  ((girth * 0.3937)^2 * (length * 0.3937)) / 800
}

locDat <- setDat %>%
  select(event, date, lat = lat2, long = long2)
indDat <- chinDat %>%
  select(fish, event, fl, circ, clip, acoustic) %>%
  mutate(size = case_when(
    fl > 75 ~ "large",
    TRUE ~ "medium"
  ),
  wt = calcWt(fl, circ)) %>% 
  left_join(., locDat, by = "event")

indDat %>% 
  group_by(clip, size) %>% 
  summarise(number = n())

# Size by clip
ggplot(indDat, aes(x = clip, y = fl)) +
  geom_boxplot() +
  samSim::theme_sleekX()

ggplot(indDat, aes(x = wt)) +
  geom_histogram(data = indDat %>% filter(clip == "Y"), fill = "red", 
                 alpha = 0.2, bins = 25) +
  geom_histogram(data = indDat %>% filter(clip == "N"), fill = "blue", 
                 alpha = 0.2, bins = 25) +
  samSim::theme_sleekX()

ggplot(indDat %>% filter(!is.na(acoustic)), 
       aes(x = fl)) +
  geom_histogram(bins = 25)
