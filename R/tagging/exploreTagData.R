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
  select(event, date, lat, long)

indDat <- chinDat %>%
  mutate(pit = paste(pitLeader, pitID, sep = ""),
         size = case_when(
           fl > 75 ~ "large",
           fl < 55 ~ "sublegal",
           TRUE ~ "medium"
         ),
         wt = calcWt(fl, circ)) %>%
  select(fish, event, fl, circ, wt, size, clip, dna, acoustic, pit) %>% 
  left_join(., locDat, by = "event")
indDat[which(indDat$pit == "NANA"), ]$pit <- NA

# write.csv(indDat, here::here("data", "taggingData", "cleanTagData.csv"),
#           row.names = FALSE)

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
