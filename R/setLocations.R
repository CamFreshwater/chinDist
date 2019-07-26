## Initial exploration of High Seas adult salmon catches
# July 22, 2019

library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
# library(PBSmapping)
# library(maps)
library(ggmap)
# library(mapdata)
library(measurements)
library(rgdal)

setDat <- read.csv(here::here("data", "taggingData", "setData.csv")) 

# Convert to decimal degrees
setDat <- setDat %>% 
  mutate(lat = paste(latDegreeStart, latMinuteStart, sep = " "),
         long = paste(longDegreeStart, longMinuteStart, sep = " ")) %>% 
  mutate(lat2 = as.numeric(measurements::conv_unit(lat, from = 'deg_dec_min', 
                                                   to = 'dec_deg')),
         long2 = -1 * as.numeric(measurements::conv_unit(long, 
                                                         from = 'deg_dec_min', 
                                                         to = 'dec_deg')))

# nAm <- map_data("world2") %>% 
#   filter(region %in% c("Canada", "USA")) %>% 
#   mutate(lat = -lat)

nAm <- ne_countries(scale = "large")

#critical habitat shape file
ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
             layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
ch84 <- spTransform(ch, CRS("+proj=longlat +datum=WGS84"))

ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  # coord_fixed(xlim = c(-128, -125), ylim = c(47, 50), ratio = 1.3) + 
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
