## Initial exploration of High Seas adult salmon catches
# July 22, 2019

library(tidyverse)
library(ggplot2)
library(PBSmapping)
library(maps)
library(ggmap)
library(mapdata)
library(measurements)

setDat <- read.csv(here::here("data", "taggingData", "tagChinData.csv")) 

# Convert to decimal degrees
setDat <- setDat %>% 
  mutate(lat = paste(latDegree, latMinute, sep = " "),
         long = paste(longDegree, longMinute, sep = " ")) %>% 
  mutate(lat2 = as.numeric(measurements::conv_unit(lat, from = 'deg_dec_min', 
                                                   to = 'dec_deg')),
         long2 = -1 * as.numeric(measurements::conv_unit(long, 
                                                         from = 'deg_dec_min', 
                                                         to = 'dec_deg')))

nAm <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))

ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-128, -125), ylim = c(47, 50), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_point(data = setDat, aes(x = long2, y = lat2),
              inherit.aes = FALSE, width = 0.1, height = 0.1) +
  geom_polygon(data = ch, aes(x = long, y = lat, group = group), 
               colour = "red", fill = NA) +
  labs(x = "Longitude", y = "Latitude") +
  samSim::theme_sleekX(legendSize = 0.8) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.2, 0.15))


ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
             layer = "SRKW_FCH_11_14_06_NAD83")
summary(ch@data)

ggplot() + 
  geom_polygon(data = ch, aes(x = long, y = lat, group = group), 
               colour = "red", fill = NA)
