## Catch data from Hargreaves sampling
# August 7, 2019

library(tidyverse)
library(ggplot2)
library(ggmap)

bark <- read.csv(here::here("data", "hargreavesData", "barkleySetLogs.csv"), 
                    stringsAsFactors = FALSE)

temp <- strsplit(bark$DATE_TIME_ISO8601, "T")

sampTime <- lapply(temp, function(x) x[2]) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rename(time = V1)

sampDate <- sapply(temp, function(x) strsplit(x[1], "-")) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rename(year = V1, month = V2, day = V3) %>% 
  cbind(., sampTime)

barkCK <- bark %>% 
  cbind(sampDate) %>% 
  select(ID_NUM:LONGITUDE, year:time, count = S_CHINOOK1, 
         countLength = S_CHINOOK2,
         countExamined = S_CHINOOK3, countBite = S_CHINOOK4,
         countClipped = S_CHINOOK5, countMarked = S_CHINOOK6) %>%
  mutate(GEAR = fct_recode(as.factor(GEAR), 
                           "PURSE SEINE" = "PURSE SEINE (TAHLOK)", 
                           "PURSE SEINE" = "PURSE SEINE (KETA)",
                           "PURSE SEINE" = "PURSE SEINE (KYNOC)",
                           "PURSE SEINE" = "PURSE SEINE (VIKING SPIRIT)"))

# How many events with juvenile chinook captured?
barkCKTrim <- barkCK %>% 
  filter(count > 0,
         !is.na(LONGITUDE),
         !is.na(LATITUDE))

barkCKTrim %>% 
  tally()

# map specs
nAm <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"), 
         long > -126, long < -124.5,
         lat > 48.75, lat < 49.5)
nAm <- ne_countries(scale = "large")
longRange = c(-126, -124.5)#range(bark$LONGITUDE, na.rm = T)
latRange = c(48.75, 49.5)#range(bark$LATITUDE, na.rm = T)

ggplot(barkCKTrim) +
  # stat_density2d(aes(x = LONGITUDE, y = LATITUDE, fill = ..level..), 
  #                geom = "polygon",
  #                alpha = 0.5) +
  geom_point(aes(x = LONGITUDE, y = LATITUDE), 
              inherit.aes = FALSE, width = 0.1, height = 0.1) + 
  # scale_fill_gradient(low = "green", high = "red") +
  # scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
  # lims(x = longRange, y = latRange) +
  geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
               color = "black", fill = "gray80") + 
  samSim::theme_sleekX()


ggplot() +
  geom_map(data=shape_map, map=shape_map, aes(long, lat, map_id=region)) +
  geom_map(
    data=filter(prop.test, season=="DJF"),
    map=shape_map, aes(fill=prop.mega, map_id=megaregion)
  )