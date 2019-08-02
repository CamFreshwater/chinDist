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
         Month %in% c("MAY", "JUN", "JUL", "AUG", "SEP"),
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
  group_by(REGION_1, HSS_REGION_1, Year) %>% 
  summarize(catch = length(REGION_1))


plotTheme <- theme(legend.title = element_text(colour = "grey30", 
                                               size = rel(0.8)),
                   panel.border = element_rect(fill = NA, colour = "grey70", 
                                               size = 1),
                   legend.key.size = unit(1, "lines"),
                   legend.text = element_text(size = rel(0.6), colour = "grey30"))

png(here::here("figs", "highSeasCatch", "catchByStock.png"), height = 2.5,
    width = 6, units = "in", res = 300)
ggplot(adDatTrim, aes(x = Year, y = catch, fill = HSS_REGION_1)) +
  geom_bar(stat = "identity") +
  plotTheme
dev.off()


# Maps by stock
adDatTrim2 <- adDat %>% 
  group_by(REGION_1) %>%
  filter(n() > 4) %>% 
  mutate(HSS_REGION_1 = fct_recode(HSS_REGION_1, 
                               OREGON_NC = "OREGON NORTH & CENTRAL")) %>% 
  ungroup()

png(here::here("figs", "highSeasCatch", "catchByStockMap.png"), height = 4,
    width = 6, units = "in", res = 300)
ggplot(data = nAm, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-129, -123), ylim = c(48, 51), ratio = 1.3) + 
  geom_polygon(color = "black", fill = "gray80") +
  geom_jitter(data = adDatTrim2, 
              aes(x = START_LONG, y = START_LAT, colour = as.factor(Year)), 
              inherit.aes = FALSE, width = 0.05, height = 0.05) + 
  labs(x = "Longitude", y = "Latitude") +
  samSim::theme_sleekX(legendSize = 0.7) +
  facet_wrap(~HSS_REGION_1)
dev.off()
