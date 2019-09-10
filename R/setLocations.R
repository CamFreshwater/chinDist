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

# Convert to decimal degrees
setDat <- setDat %>% 
  mutate(lat = paste(latDegreeStart, latMinuteStart, sep = " "),
         long = paste(longDegreeStart, longMinuteStart, sep = " ")) %>% 
  mutate(lat2 = as.numeric(measurements::conv_unit(lat, from = 'deg_dec_min', 
                                                   to = 'dec_deg')),
         long2 = -1 * as.numeric(measurements::conv_unit(long, 
                                                         from = 'deg_dec_min', 
                                                         to = 'dec_deg')),
         julDay = as.POSIXlt(date, format = "%d/%m/%Y")$yday) %>% 
  select(event, set, date, julDay, lat = lat2, long = long2, time)

# write.csv(setDat, here::here("data", "taggingData", "cleanSetData.csv"),
#           row.names = FALSE)

nAm <- ne_countries(scale = "large")

#critical habitat shape file
ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
             layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
ch84 <- spTransform(ch, CRS("+proj=longlat +datum=WGS84"))


## Plot set locations
setMap <- ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-127, -124.5), ylim = c(48, 49.25), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_point(data = setDat, aes(x = long, y = lat, fill = julDay),
             inherit.aes = FALSE, shape = 21) +
  viridis::scale_fill_viridis(option = "magma") +
  geom_polygon(data = ch84, aes(x = long, y = lat, group = group), 
               colour = "red", fill = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Julian Day") +
  samSim::theme_sleekX(legendSize = 0.8) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.1, 0.2))
# saveRDS(setMap, here::here("generatedData", "setMap.RDS"))


png(here::here("figs", "maps", "setMap.png"), height = 5, width = 7,
    units = "in", res = 400)
setMap
dev.off()
