## Examine juvenile catches from IPES survey to identify potential receiver 
# locations for NOAA
# Aug. 1, 2019

library(tidyverse)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggmap)
library(measurements)
library(rgdal)

ipesDat <- read.csv(here::here("data", "highSeas", "ipesChinook2017_19.csv")) %>% 
  rename(lat = START_LATITUDE, long = START_LONGITUDE, year = TRIP_YEAR)

## Plot set locations to constrain subsequent analyses
nAm <- ne_countries(scale = "large")

# trim to near Brooks
trimDat <- ipesDat %>% 
  filter(lat > 49.75, lat < 50.25,
         long < -127.3, long > -128.25) 
  

p <- ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-128.75, -123), ylim = c(48.25, 52), ratio = 1.3) +
  geom_polygon(color = "black", fill = "gray80") +
  labs(x = "Longitude", y = "Latitude") +
  samSim::theme_sleekX(legendSize = 0.8) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.2, 0.15))

p + 
  geom_point(data = trimDat, aes(x = long, y = lat, color = as.factor(year)),
             inherit.aes = FALSE)
  
summDat <- ipesDat %>% 
  group_by(EVENT_DATE) %>% 
  summarize(catch = n()) %>% 
  mutate(highCatch = case_when(
    catch > 4 ~ "yes",
    TRUE ~ "no"
  )) 

catchDat <- ipesDat %>% 
  select(EVENT_DATE, lat, long) %>% 
  distinct() %>% 
  left_join(., summDat, by = "EVENT_DATE")
  
ggplot(catchDat) +
  stat_density2d(aes(x = long, y = lat, fill = ..level..), 
                 geom = "polygon",
                 alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  geom_point(aes(x = long, y = lat)) +
  lims(x = c(-129.5, -123), y = c(48, 52)) +
  geom_polygon(data = nAm, aes(x = long, y = lat, group = group), 
               color = "black", fill = "gray80") 
## Evidence of high catches near Lippey Point