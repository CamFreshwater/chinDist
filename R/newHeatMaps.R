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

## Calculate CPUE for each set, accounting for clip vs. unclipped
blankCount1 <- setDat %>%
  select(event) %>%
  mutate(clip = "Y")
blankCount <- blankCount1 %>%
  mutate(clip = "N") %>%
  rbind(., blankCount1)

catch <- chinDat %>%
  group_by(event, clip) %>%
  summarise(nChin = n()) %>%
  full_join(., blankCount, by = c("event", "clip")) %>%
  replace_na(list(nChin = 0)) %>% 
  full_join(setDat, ., by = "event") %>% 
  arrange(event) %>%
  select(set, lat, long, clip, catch = nChin)

# convert catch data to spatial
setPts <- catch %>% 
  group_by(set, lat, long) %>% 
  summarize(totalCatch = sum(catch)) %>% 
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%  
  # st_transform(., crs = 32610) %>%
  #bind coordinates as numeric
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
# %>%
#   #clip to habitat shapefile
#   st_intersection(., hab)

hex_polygons$effort <- lengths(st_intersects(hex_polygons, setPts))

zero <- hex_polygons %>% 
  filter(effort == 0)
nonZero <- hex_polygons %>% 
  filter(effort > 0)

p2 <- ggplot() +
  geom_sf(data = hex_polygons) +
  # geom_sf(data = zero, aes(fill = "white")) %>% 
  geom_sf(data = nonZero, aes(fill = effort)) +
  scale_fill_viridis_c(option = "magma", "# of Sets") +
  theme_void() 

plot(p2)
