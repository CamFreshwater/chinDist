## Initial exploration of tagging data
# July 22, 2019

library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggmap)
library(measurements)
library(rgdal)

setDat <- read.csv(here::here("data", "taggingData", "setData.csv")) 
chinDat <- read.csv(here::here("data", "taggingData", "tagChinData.csv"))

# Convert to decimal degrees
setDat <- setDat %>% 
  mutate(lat = paste(latDegreeStart, latMinuteStart, sep = " "),
         long = paste(longDegreeStart, longMinuteStart, sep = " ")) %>% 
  mutate(lat2 = as.numeric(measurements::conv_unit(lat, from = 'deg_dec_min', 
                                                   to = 'dec_deg')),
         long2 = -1 * as.numeric(measurements::conv_unit(long, 
                                                         from = 'deg_dec_min', 
                                                         to = 'dec_deg')))

nAm <- ne_countries(scale = "large")

#critical habitat shape file
ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
             layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
ch84 <- spTransform(ch, CRS("+proj=longlat +datum=WGS84"))


## Plot set locations
setMap <- ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-128.5, -123), ylim = c(48.25, 51), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_point(data = setDat, aes(x = long2, y = lat2),
             inherit.aes = FALSE, width = 0.1, height = 0.1) +
  geom_polygon(data = ch84, aes(x = long, y = lat, group = group), 
               colour = "red", fill = NA) +
  labs(x = "Longitude", y = "Latitude") +
  samSim::theme_sleekX(legendSize = 0.8) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.2, 0.15))

# png(here::here("figs", "maps", "setMap.png"), height = 5, width = 7, 
#     units = "in", res = 400)
setMap
# dev.off()

#-----

## Plot heat map of catch data

# Sum chinook catch data
count <- chinDat %>%
  group_by(event) %>%
  summarise(nChin = n())
# Add lat/longs from set data to catch data
catchDat <- setDat %>%
  select(event, date, lat = lat2, long = long2) %>%
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
               color = "black", fill = "gray80") +
  geom_polygon(data = ch84, aes(x = long, y = lat, group = group), 
               colour = "red", fill = NA)

# ALT VERSION (doesn't account for effort as well)
# Merge catch and fish data
locDat <- setDat %>%
  select(event, date, lat = lat2, long = long2)
indDat <- chinDat %>%
  select(fish, event, fl, circ, clip) %>% 
  mutate(size = case_when(
    fl > 75 ~ "large",
    TRUE ~ "medium"
  ))  %>% 
  left_join(., locDat, by = "event")

p <- ggplot(indDat) +
  stat_density_2d(aes(x = long, y = lat, fill = ..level..),
                  geom = "polygon",
                  alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  geom_point(data = locDat, aes(x = long, y = lat)) +
  lims(x = c(-126.3, -125.3), y = c(48.35, 49.25)) +
  coord_fixed(ratio = 1.3) +
  geom_polygon(data = nAm, aes(x = long, y = lat, group = group), 
               color = "black", fill = "gray80") +
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
