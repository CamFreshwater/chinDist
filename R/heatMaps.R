## Heat maps of raw catch data
#Sep. 10, 2019

library(tidyverse)
library(ggplot2)
library(maps)

## Sampling data
setDat <- read.csv(here::here("data", "taggingData", "cleanSetData.csv")) 
chinDat <- read.csv(here::here("data", "taggingData", "cleanTagData.csv"))


wCan <- map_data("world", region = "canada") %>%
  filter(long < -110)

chinHeat <- ggplot(chinDat) +
  stat_density2d(aes(x = long, y = lat, fill = stat(level)), 
                 geom = "polygon",
                 alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  lims(x = c(-130, -124), y = c(48, 52)) +
  geom_polygon(data = wCan, aes(x = long, y = lat, group = group),
               color = "black", fill = "gray80") +
  geom_point(data = setDat, aes(x = long, y = lat), alpha = 0.4) +
  coord_fixed(xlim = c(-126.3, -125.3), ylim = c(48.35, 49.2), ratio = 1.3,
              expand = TRUE) +
  labs(x = "Longitude", ylab = "Latitude") +
  samSim::theme_sleekX()

chinHeat +
  facet_wrap(~size)

## First pass on version that better accounts for zero sets, but need to double 
# check how effort is rolled in

count <- chinDat %>%
  group_by(event, clip) %>%
  summarise(nChin = n())
# Add lat/longs from set data to catch data
catchDat <- setDat %>%
  select(event, date, lat, long) %>%
  left_join(., count, by = "event") %>%
  replace_na(list(nChin = 0))

chinHeat2 <- ggplot(catchDat, aes(x = long, y = lat, color = nChin)) +
  geom_point(aes(color = nChin), shape = "") +
  stat_density2d(aes(fill = ..level..), 
                 # n = 100, 
                 contour = TRUE, 
                 geom = "polygon") +
  scale_color_continuous(guide = FALSE) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  lims(x = c(-130, -124), y = c(48, 52)) +
  geom_polygon(data = wCan, aes(x = long, y = lat, group = group),
               color = "black", fill = "gray80") +
  # geom_point(data = setDat, aes(x = long, y = lat), inherit.aes = FALSE) +
  coord_fixed(xlim = c(-126.3, -125.3), ylim = c(48.35, 49.2), ratio = 1.3,
              expand = TRUE) +
  labs(x = "Longitude", ylab = "Latitude") +
  samSim::theme_sleekX()

## Old buggy version with polygon issues

## Spatial files
# nAm <- ne_countries(scale = "large")
# #critical habitat shape file
# ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
#               layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
# ch84 <- spTransform(ch, CRS("+proj=longlat +datum=WGS84"))

## Plot heat map of catch data

# Sum chinook catch data
# count <- chinDat %>%
#   group_by(event) %>%
#   summarise(nChin = n())
# # Add lat/longs from set data to catch data
# catchDat <- setDat %>%
#   select(event, date, lat, long) %>%
#   left_join(., count, by = "event") %>%
#   replace_na(list(nChin = 0))
# 
# ggplot(catchDat) +
#   stat_density2d(aes(x = long, y = lat, fill = stat(level)), 
#                  geom = "polygon",
#                  alpha = 0.5) +
#   scale_fill_gradient(low = "green", high = "red") +
#   scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
#   geom_point(aes(x = long, y = lat)) +
#   geom_polygon(data = nAm, aes(x = long, y = lat, group = group), 
#                color = "black", fill = "gray80") +
#   lims(x = c(-126.3, -125.3), y = c(48.35, 49.4))
# 
# 
# # ALT VERSION (doesn't account for effort as well)
# # Merge catch and fish data
# ggplot(chinDat) +
#   stat_density_2d(aes(x = long, y = lat, fill = ..level..),
#                   geom = "polygon",
#                   alpha = 0.5) +
#   scale_fill_gradient(low = "green", high = "red") +
#   scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
#   geom_point(data = setDat, aes(x = long, y = lat)) +
#   coord_map(xlim = c(-126.3, -125.3), ylim = c(48, 49.5)) +
#   geom_polygon(data = nAm, aes(x = long, y = lat, group = group), 
#              color = "black", fill = "gray80", inherit.aes = FALSE)
#   coord_cartesian(xlim = c(-126.3, -125.3), ylim = c(48.35, 49.5))
#   # lims(x = c(-126.3, -125.3), y = c(48.35, 49.5)) 
