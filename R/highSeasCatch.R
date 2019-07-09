## Initial exploration of High Seas adult salmon catches
# July 5, 2019

library(tidyverse)
library(ggplot2)
library(PBSmapping)
library(maps)
library(ggmap)
library(mapdata)


# Restrict to stock ID'd samples
adDatFull <- read.csv(here::here("data", "highSeas", "adultChinookHS.csv")) 

str_detect(adDatFull$STATION, "GS")

adDat <- adDatFull %>% 
  filter(PROB_1 > 0.5, 
         # remove non-WCVI capture locations
         START_LAT < 51,
         !str_detect(STATION, "BR"),
         !str_detect(STATION, "JS"),
         !str_detect(STATION, "QCS"),
         !str_detect(STATION, "PS"),
         !str_detect(STATION, "GI"),
         !str_detect(STATION, "GS"),
         !str_detect(STATION, "HW")
         )
unique(adDat$STATION)

psDat <- adDat %>% 
  filter(REGION_1 == "PUGET SOUND")

# Plot capture locations
nAm <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))

ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-138, -119), ylim = c(48, 59.5), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_jitter(data = adDat, aes(x = START_LONG, y = START_LAT), 
              inherit.aes = FALSE, width = 0.1, height = 0.1) + 
  labs(x = "Longitude", y = "Latitude") +
  samSim::theme_sleekX(legendSize = 0.8) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(0.2, 0.15))

# Stacked bar plots
adDatTrim <- adDat %>% 
  group_by(REGION_1) %>%
  filter(n() > 4) %>% 
  group_by(REGION_1, Year) %>% 
  summarize(catch = length(REGION_1))

ggplot(adDatTrim, aes(x = Year, y = catch, fill = REGION_1)) +
  geom_bar(stat = "identity")
