## Heat maps of raw catch data
#Sep. 10, 2019

library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggmap)
library(measurements)
library(rgdal)

## Sampling data
setDat <- read.csv(here::here("data", "taggingData", "cleanSetData.csv")) 
chinDat <- read.csv(here::here("data", "taggingData", "cleanTagData.csv"))

## Spatial files
nAm <- ne_countries(scale = "large")
#critical habitat shape file
ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
              layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
ch84 <- spTransform(ch, CRS("+proj=longlat +datum=WGS84"))


## Plot heat map of catch data

# Sum chinook catch data
count <- chinDat %>%
  group_by(event) %>%
  summarise(nChin = n())
# Add lat/longs from set data to catch data
catchDat <- setDat %>%
  select(event, date, lat, long) %>%
  left_join(., count, by = "event") %>%
  replace_na(list(nChin = 0))

ggplot(catchDat) +
  stat_density2d(aes(x = long, y = lat, fill = ..level..), 
                 geom = "polygon",
                 alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  geom_point(aes(x = long, y = lat)) +
  lims(x = c(-126.3, -125.3), y = c(48.35, 49.2)) +
  geom_polygon(data = nAm, aes(x = long, y = lat, group = group), 
               color = "black", fill = "gray80") 

# ALT VERSION (doesn't account for effort as well)
# Merge catch and fish data
p <- ggplot(chinDat) +
  stat_density_2d(aes(x = long, y = lat, fill = ..level..),
                  geom = "polygon",
                  alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  geom_point(data = setDat, aes(x = long, y = lat)) +
  geom_polygon(data = nAm, aes(x = long, y = lat, group = group), 
               color = "black", fill = "gray80") +
  # lims(x = c(-126.3, -125.3), y = c(48.35, 49.5)) +
  lims(x = c(-126.5, -124.5), y = c(48, 50)) + 
  coord_fixed(ratio = 1.3) +
  samSim::theme_sleekX()


ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-127, -124.5), ylim = c(48, 49.25), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  # geom_point(data = setDat, aes(x = long, y = lat),
  #            inherit.aes = FALSE, shape = 21) +
  stat_density_2d(data = chinDat, aes(x = long, y = lat, fill = ..level..),
                  geom = "polygon",
                  alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  labs(x = "Longitude", y = "Latitude", fill = "Julian Day") +
  samSim::theme_sleekX() 





png(here::here("figs", "maps", "hatcheryMap.png"), height = 5, width = 7,
    units = "in", res = 400)
p +
  facet_wrap(~clip)
dev.off()


png(here::here("figs", "maps", "sizeMap.png"), height = 5, width = 7,
    units = "in", res = 400)
p +
  facet_wrap(~size)
dev.off()
