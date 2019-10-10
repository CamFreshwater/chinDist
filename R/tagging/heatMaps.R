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

# convert set data to spatial (this crs is in m, can also use 4326 for WGS84 
# projection)
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
hex_points <- spsample(as_Spatial(hab), type = "hexagonal", cellsize = 0.05)
hex_polygons <- HexPoints2SpatialPolygons(hex_points) %>%
  st_as_sf(crs = st_crs(setPts))

# check dimensions
# temp_points <- spsample(as_Spatial(hab), type = "regular", cellsize = 0.04)
# coordinates(temp_points)[1:10, ]

# NOTE eventually add clip vs. unclip, size classes and hatcheries
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

# Incorporate clustering (basically accounts for spatial autocorrelation 
# associated with a set integrating over fairly large distances)

library(spdep)

#find the centroid of each hexagon and convert to a matrix of points
hex_points <- do.call(rbind, st_geometry(st_centroid(hex_polygons))) %>%
  unlist() %>%
  as.matrix.data.frame()

#use a k-nearest-neighbour algorithm to find which shape neighbour which
#super easy for the hexagons analytically obvs but important for e.g. using the 
#arbitrary discreat boundaries instead
neighbouring_hexes <- knn2nb(knearneigh(hex_points, k = 6), 
                             row.names = rownames(hex_points)) %>%
  include.self()

#calculate the local G for a given variable (#bus stops) using the neihbours found previously
localGvalues <- localG(x = as.numeric(hex_polygons$catch),
                       listw = nb2listw(neighbouring_hexes, style = "B"),
                       zero.policy = TRUE)

#bind this back to the sf as a numeric variable column
hex_polygons <- hex_polygons %>% 
  mutate(smoothed_catch = as.numeric(localGvalues),
         smoothed_cpue = smoothed_catch / effort)

nonZeroSmooth <- hex_polygons %>% 
  filter(effort > 0) 

#plot the statistic
#+ve indicates more than would be expected
p3 <- ggplot() +
  geom_sf(data = hex_polygons) +
  geom_sf(data = nonZeroSmooth, aes(fill = smoothed_cpue)) +
  scale_fill_viridis_c(option = "magma", "Smoothed\nCPUE") +
  theme_void() 
plot(p3)
