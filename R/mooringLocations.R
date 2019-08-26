## Mooring locations
# Aug 26, 2019

library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggmap)
library(measurements)
library(rgdal)
library(viridis)

nAm <- ne_countries(scale = "large")

#critical habitat shape file
ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
              layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
ch84 <- spTransform(ch, CRS("+proj=longlat +datum=WGS84"))


## Plot set locations
setMap <- ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-128, -123.5), ylim = c(48, 51), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  # geom_point(data = setDat, aes(x = long2, y = lat2, fill = julDay),
  #            inherit.aes = FALSE, shape = 21) +
  # scale_fill_viridis(option = "magma") +
  geom_polygon(data = ch84, aes(x = long, y = lat, group = group), 
               colour = "red", fill = NA) +
  labs(x = "Longitude", y = "Latitude") +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.1, 0.2))
