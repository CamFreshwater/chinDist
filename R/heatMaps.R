## Heat maps of catch per unit effort data
#Sep. 20, 2019

#Second pass using tips from 
#http://www.robert-hickman.eu/post/getis-ord-heatmaps-tutorial/

library(ggplot2)
library(ggmap)
library(measurements)
library(rgdal)
library(sf)
library(raster)
library(tidyverse)

## First calculate catch for set data
setDat <- read.csv(here::here("data", "taggingData", "cleanSetData.csv")) 
chinDat <- read.csv(here::here("data", "taggingData", "cleanTagData.csv"))

# convert set data to spatial
setPts <- setDat %>% 
  select(event, lat, long) %>% 
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%  
  cbind(., st_coordinates(.))

#convert ind. data to spatial
chinTrim <- chinDat %>%
  select(event, lat, long, clip) %>% 
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%  
  cbind(., st_coordinates(.))

# habitat shape file
hab <- st_read(
  here::here("data", "criticalHabitatShapeFiles",
             "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016.shp")) %>% 
  st_transform(4326)

# create equidistant bins from shape file
hex_points <- spsample(as_Spatial(hab), type = "hexagonal", cellsize = 0.04)
hex_polygons <- HexPoints2SpatialPolygons(hex_points) %>%
  st_as_sf(crs = st_crs(setPts))

# eventually add clip vs. unclip, size classes and hatcheries
hex_polygons <- hex_polygons %>% 
  mutate(effort = lengths(st_intersects(hex_polygons, setPts)),
         catch = lengths(st_intersects(hex_polygons, chinTrim)),
         cpue = catch / effort)

nonZero <- hex_polygons %>% 
  filter(effort > 0) 

p2 <- ggplot() +
  geom_sf(data = hex_polygons) +
  geom_sf(data = nonZero, aes(fill = cpue)) +
  scale_fill_viridis_c(option = "magma", "CPUE") +
  theme_void() 
plot(p2)

# ggplot() +
#   geom_sf(data = hex_polygons) +
#   geom_sf(data = nonZero, aes(fill = effort)) +
#   scale_fill_viridis_c(option = "magma", "effort") +
#   theme_void()

