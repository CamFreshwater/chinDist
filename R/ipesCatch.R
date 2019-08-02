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
nAm <- ne_countries(scale = "medium")

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
  
summDat <- trimDat %>% 
  group_by(EVENT_DATE) %>% 
  summarize(catch = n()) %>% 
  mutate(highCatch = case_when(
    catch > 4 ~ "yes",
    TRUE ~ "no"
  ))

trimDat2 <- left_join(trimDat, summDat, by = "EVENT_DATE")

p + 
  geom_point(data = trimDat2, aes(x = long, y = lat, 
                                  color = as.factor(highCatch)),
             inherit.aes = FALSE)

## Evidence of high catches near Lippey Point