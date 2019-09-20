## Initial exploration of tagging data
# July 22, 2019

library(tidyverse)
library(ggplot2)
library(ggmap)
library(measurements)
library(rgdal)
library(sf)
library(raster)

setDat1 <- read.csv(here::here("data", "taggingData", "setData.csv")) 

# Convert to decimal degrees
setDatFull <- setDat1 %>% 
  mutate(startLat = paste(latDegreeStart, latMinuteStart, sep = " "),
         startLong = paste(longDegreeStart, longMinuteStart, sep = " ")) %>% 
  mutate(startLat = as.numeric(conv_unit(startLat, from = 'deg_dec_min',
                                         to = 'dec_deg')),
         startLong = -1 * as.numeric(conv_unit(startLong, from = 'deg_dec_min',
                                               to = 'dec_deg')),
         julDay = as.POSIXlt(date, format = "%d/%m/%Y")$yday)

wCan <- map_data("world", region = "canada") %>%
  filter(long < -110)

#critical habitat shape file
ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
             layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
#transform to WGS84
ch84 <- sp::spTransform(ch, CRS("+proj=longlat +datum=WGS84"))

## Plot set locations
setMap <- ggplot(data = wCan, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-127, -124.5), ylim = c(48, 49.25), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_point(data = setDatFull, aes(x = startLong, y = startLat, fill = julDay),
             inherit.aes = FALSE, shape = 21) +
  viridis::scale_fill_viridis(option = "magma") +
  geom_polygon(data = ch84, aes(x = long, y = lat, group = group), 
               colour = "red", fill = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Julian Day") +
  samSim::theme_sleekX(legendSize = 0.8) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.1, 0.2))
setMap
# saveRDS(setMap, here::here("generatedData", "setMap.RDS"))

png(here::here("figs", "maps", "setMap.png"), height = 5, width = 7,
    units = "in", res = 400)
setMap
dev.off()


### Generate grid based on average set distance
dist_geo <- function(lat_a, lon_a, lat_b, lon_b) { 
  library(geosphere)
  distm(c(lon_a, lat_a), c(lon_b, lat_b), fun = distGeo) 
} 

distDat <- setDatFull %>% 
  filter(!is.na(latDegreeEnd)) %>% 
  mutate(endLat = paste(latDegreeEnd, latMinuteEnd, sep = " "),
         endLong = paste(longDegreeEnd, longMinuteEnd, sep = " ")) %>% 
  mutate(endLat = as.numeric(conv_unit(endLat, from = 'deg_dec_min',
                                       to = 'dec_deg')),
         endLong = -1 * as.numeric(conv_unit(endLong, from = 'deg_dec_min',
                                             to = 'dec_deg')),
         setDist = mapply(lat_a = startLat, lon_a = startLong,
                          lat_b = endLat, lon_b = endLong,
                          FUN = dist_geo)) %>% 
  dplyr::select(set, endLat, endLong, setDist)

distDat %>%
  filter(!setDist > 5000) %>%
  summarize(mean = mean(setDist),
            sd = sd(setDist),
            max = max(setDist),
            min = min(setDist))
#start w/ 5

hab <- st_read(
  here::here("data", "criticalHabitatShapeFiles",
             "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016.shp")) %>% 
  st_transform(32610)

setPts <- setDatFull %>% 
  dplyr::select(set, startLat, startLong) %>% 
  st_as_sf(., coords = c("startLong", "startLat"), crs = 4326) %>%  
  st_transform(., crs = 32610)  
#set CRS to 4326 to match ch84 file which is WGS84, then convert to UTM for 
#finer scale res

#make grid
grid75 <- st_make_grid(hab, cellsize = c(7500, 7500)) %>% 
  st_sf(gridID = 1:length(.))

# make labels
grid_lab <- st_centroid(grid75) %>% 
  cbind(st_coordinates(.))

griddedMap <- ggplot() +
  geom_sf(data = hab, fill = 'white', lwd = 0.05) +
  geom_sf(data = setPts, color = 'red', size = 1.7) +
  geom_sf(data = grid75, fill = 'transparent', lwd = 0.3) +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = gridID), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "")

# png(here::here("figs", "maps", "griddedMap_7.5Cells.png"), height = 5, 
#     width = 7, units = "in", res = 400)
griddedMap
# dev.off()

# id which grid different points are in
gridCells <- setPts %>% 
  st_join(grid75, join = st_intersects) %>% 
  as.data.frame()

setDat <- setDatFull %>%
  left_join(., distDat, by = "set") %>% 
  left_join(., gridCells, by = "set") %>% 
  dplyr::select(event, set, time, date, julDay, gridID, lat = startLat, 
                long = startLong, setDist)

# write.csv(setDat, here::here("data", "taggingData", "cleanSetData.csv"),
#           row.names = FALSE)
