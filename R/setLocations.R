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

setDat1 <- read.csv(here::here("data", "taggingData", "setData.csv")) 

# Convert to decimal degrees
setDatFull <- setDat1 %>% 
  mutate(startLat = paste(latDegreeStart, latMinuteStart, sep = " "),
         startLong = paste(longDegreeStart, longMinuteStart, sep = " ")) %>% 
  mutate(startLat = as.numeric(measurements::conv_unit(startLat, 
                                                       from = 'deg_dec_min',
                                                       to = 'dec_deg')),
         startLong = -1 * as.numeric(measurements::conv_unit(startLong, 
                                                         from = 'deg_dec_min', 
                                                         to = 'dec_deg')),
         julDay = as.POSIXlt(date, format = "%d/%m/%Y")$yday)

setDat <- setDatFull %>% 
  select(event, set, date, julDay, lat = startLat, long = startLong, time)

# write.csv(setDat, here::here("data", "taggingData", "cleanSetData.csv"),
#           row.names = FALSE)

nAm <- ne_countries(scale = "large")

#critical habitat shape file
ch <- readOGR(here::here("data", "criticalHabitatShapeFiles"), 
             layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016")
#transform to WGS84
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


### Generate grid based on average set distance
dist_geo <- function(lat_a, lon_a, lat_b, lon_b) { 
  library(geosphere)
  distm(c(lon_a, lat_a), c(lon_b, lat_b), fun = distGeo) 
} 

distDat <- setDatFull %>% 
  filter(!is.na(latDegreeEnd)) %>% 
  mutate(endLat = paste(latDegreeEnd, latMinuteEnd, sep = " "),
         endLong = paste(longDegreeEnd, longMinuteEnd, sep = " ")) %>% 
  mutate(endLat = as.numeric(measurements::conv_unit(endLat,
                                                from = 'deg_dec_min',
                                                to = 'dec_deg')),
         endLong = -1 * as.numeric(measurements::conv_unit(endLong,
                                                    from = 'deg_dec_min',
                                                    to = 'dec_deg')),
         setDist = mapply(lat_a = startLat, lon_a = startLong,
                          lat_b = endLat, lon_b = endLong,
                          FUN = dist_geo)) %>% 
  select(event, set, date, time, julDay, startLat, startLong, endLat, endLong,
         setDist)

# distDat %>% 
#   filter(!setDist > 5000) %>% 
#   summarize(mean = mean(setDist),
#             sd = sd(setDist),
#             max = max(setDist),
#             min = min(setDist))
#start w/ 5

wCan <- map_data("world", region = "canada") %>%
  filter(long < -110)

# sets <- 
ggplot(distDat) +
  geom_point(aes(x = startLong, y = startLat), alpha = 0.4) +
  lims(x = c(-130, -124), y = c(48, 52)) +
  geom_polygon(data = wCan, aes(x = long, y = lat, group = group),
               color = "black", fill = "gray80") +
  coord_fixed(xlim = c(-126.3, -125.3), ylim = c(48.35, 49.2), ratio = 1.3,
              expand = TRUE) +
  labs(x = "Longitude", y = "Latitude", fill = "Relative\nAbundance") +
  theme_bw()

  
## Option 1:  
library(sf)
library(raster)

grid <- ch_sf %>% 
  st_make_grid(cellsize = 0.1, what = "centers") %>% # grid of points
  st_intersection(ch_sf)     

ggplot() +
  geom_sf(data = ch_sf) +
  geom_sf(data = grid) +
  geom_point(data = distDat, aes(x = startLong, y = startLat), alpha = 0.4)

## Option 2:



hab <- st_read(
  here::here("data", "criticalHabitatShapeFiles",
             "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016.shp")) %>% 
  st_transform(32610)

setPts <- distDat %>% 
  select(set, startLat, startLong) %>% 
  st_as_sf(., coords = c("startLong", "startLat"), crs = 32610) 
#set CRS to 4326 to match ch84 file which is WGS84 (shape file originally in 
#NAD83)

#make grid
grid_5 <- st_make_grid(hab, cellsize = c(500,0 5000)) %>% 
  st_sf(grid_id = 1:length(.))

# make labels
grid_lab <- st_centroid(grid_5) %>% 
  cbind(st_coordinates(.))

ggplot() +
  geom_sf(data = hab, fill = 'white', lwd = 0.05) +
  geom_sf(data = setPts, color = 'red', size = 1.7) + 
  geom_sf(data = grid_5, fill = 'transparent', lwd = 0.3) +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "")



library(sf)
library(ggplot2)

# read nc polygon data and transform to UTM 
nc <- st_read(system.file('shape/nc.shp', package = 'sf')) %>%
  st_transform(32617)

# random sample of 5 points
pts <- st_sample(nc, size = 5) %>% st_sf

# create 50km grid - here you can substitute 200 for 50000
grid_50 <- st_make_grid(nc, cellsize = c(50000, 50000)) %>% 
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
grid_lab <- st_centroid(grid_50) %>% 
  cbind(st_coordinates(.))

# view the sampled points, polygons and grid
ggplot() +
  geom_sf(data = nc, fill = 'white', lwd = 0.05) +
  geom_sf(data = pts, color = 'red', size = 1.7) + 
  geom_sf(data = grid_50, fill = 'transparent', lwd = 0.3) +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "")

# which grid square is each point in?
pts %>% st_join(grid_50, join = st_intersects) %>% as.data.frame








### read shapefile
library(rgdal)
shp <- readOGR(here::here("data", "exampleData", "nybb_13a"), layer = "nybb")

proj4string(shp)  # units us-ft

### define coordinates and convert to SpatialPointsDataFrame
poi <- data.frame(x=c(919500, 959500, 1019500, 1049500, 1029500, 989500),
                  y=c(130600, 150600, 180600, 198000, 248000, 218000),
                  id="A", stringsAsFactors=F)
coordinates(poi) <- ~ x + y
proj4string(poi) <- proj4string(shp)

### define SpatialGrid object
bb <- bbox(shp)
cs <- c(3.28084, 3.28084)*6000  # cell size 6km x 6km (for illustration)
# 1 ft = 3.28084 m
cc <- bb[, 1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
grd
# cellcentre.offset 923018 129964
# cellsize           19685  19685
# cells.dim              8      8

sp_grd <- SpatialGridDataFrame(grd,
                               data=data.frame(id=1:prod(cd)),
                               proj4string=CRS(proj4string(shp)))
summary(sp_grd)

