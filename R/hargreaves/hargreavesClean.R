## Clean catch data from Hargreaves sampling in Barkley Sound
# August 7, 2019
# See text files for details

library(tidyverse)
library(ggplot2)


### CLEAN SET DATA

bark <- read.csv(here::here("data", "hargreaves", "barkleySetLogs.csv"), 
                    stringsAsFactors = FALSE)

temp <- strsplit(bark$DATE_TIME_ISO8601, "T")

sampTime <- lapply(temp, function(x) x[2]) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rename(time = V1)

julDay <- sapply(temp, function(x) as.POSIXlt(x[1], format = "%Y-%m-%d")$yday)

sampDate <- sapply(temp, function(x) strsplit(x[1], "-")) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  mutate(julianDay = julDay) %>% 
  rename(year = V1, month = V2, day = V3) %>% 
  cbind(., sampTime)

barkSets <- bark %>% 
  cbind(sampDate) %>% 
  filter(!is.na(month)) %>% 
  select(ID_NUM:LONGITUDE, year:time, count = S_CHINOOK1, 
         countLength = S_CHINOOK2,
         countExamined = S_CHINOOK3, countBite = S_CHINOOK4,
         countClipped = S_CHINOOK5, countMarked = S_CHINOOK6) %>%
  mutate(GEAR = fct_recode(as.factor(GEAR), 
                           "PURSE SEINE" = "PURSE SEINE (TAHLOK)", 
                           "PURSE SEINE" = "PURSE SEINE (KETA)",
                           "PURSE SEINE" = "PURSE SEINE (KYNOC)",
                           "PURSE SEINE" = "PURSE SEINE (VIKING SPIRIT)"),
         month = as.factor(month),
         chinookPresent = case_when(
           count > 0 ~ "Y",
           count == 0 ~ "N"
         ))

sampDate2 <- barkSets %>% 
  select(ID_NUM, year:time)

# How many events with juvenile chinook captured?
barkSets %>% 
  filter(count > 0) %>% 
  tally()

barkSets %>% 
  group_by(GEAR) %>% 
  tally()

## Sampling breakdown
plotGear <- barkSets %>% 
  group_by(GEAR, month, chinookPresent) %>% 
  tally(name = "nn")

ggplot(plotGear, aes(x = as.factor(month), y = nn, fill = chinookPresent)) +
  geom_bar(stat = "identity") +
  facet_wrap(~GEAR)
# majority capture by purse seine and CMBS (beach seine) in May-July

barkSets %>% 
  filter(GEAR %in% c("100' CMBS", "150' CMBS", "PURSE SEINE")) %>% 
  select(PURPOSE) %>% 
  distinct()
# ABORTED, ID ERROR, REPEATED and EXPLORATORY should likely be removed

barkSetsTrim <- barkSets %>% 
  filter(GEAR %in% c("100' CMBS", "150' CMBS", "PURSE SEINE"),
         !PURPOSE %in% c("ABORTED", "ID ERROR, REPEATED"),
         !is.na(LATITUDE)) %>% 
  select(id = ID_NUM, gear = GEAR, purpose = PURPOSE, location = LOCATION,
         lat = LATITUDE, long = LONGITUDE, year:chinookPresent)

write.csv(barkSetsTrim, here::here("data", "hargreaves", 
                                 "barkleySets_cleaned.csv"), row.names = FALSE)


### CLEAN FISH LENGTH AND CWT DATA
barkFish <- read.csv(here::here("data", "hargreaves", "barkleyNonPredator.csv"),
                     stringsAsFactors = FALSE, na.strings = "")
barkCWT <- read.csv(here::here("data", "hargreaves", "barkleyCWTRecovery.csv"),
                    stringsAsFactors = FALSE, na.strings = "")

#sets kept based on gear and comments
goodSets <- unique(barkSetsTrim$id)

barkFishTrim <- barkFish %>% 
  filter(ID_NUM %in% goodSets,
         SPECIES == "CHINOOK SALMON SMOLT") %>% 
  left_join(., sampDate2, by = "ID_NUM") %>% 
  mutate(fishNumber = paste(ID_NUM, FISH_NUM, sep = "_")) %>% 
  select(id = ID_NUM, fishNumber, species = SPECIES,
         lat = LATITUDE, long = LONGITUDE, year:time, fl = LENGTH, 
         weight = WEIGHT, bitemark = BITEMARK, CWT, marked = MARKED, 
         pres = PRESERV)

write.csv(barkFishTrim, here::here("data", "hargreaves", 
                                   "barkleyFish_cleaned.csv"), 
          row.names = FALSE)


minRelJulian <- as.POSIXlt(barkCWT$FIRST_REL, format = "%m/%d/%Y")$yday
# maxRelJulian <- barkCWT$LAST_REL %>% 
  # strsplit(., "/") %>% 
  # do.call(rbind, .) %>% 
  # as.data.frame() %>% 
  # select(minRelMonth = V1, minRelDay = V2)

barkCWTTrim <- barkCWT %>% 
  mutate(fishNumber = paste(ID_NUM, FISH_NUM, sep = "_"),
         minRelJulian = as.POSIXlt(FIRST_REL, format = "%m/%d/%Y")$yday,
         maxRelJulian = as.POSIXlt(LAST_REL, format = "%m/%d/%Y")$yday) %>% 
  filter(ID_NUM %in% goodSets,
         SPECIES == "CHINOOK SALMON SMOLT") %>% 
  left_join(., sampDate2, by = "ID_NUM") %>% 
  select(id = ID_NUM, fishNumber, species = SPECIES,
         lat = LATITUDE, long = LONGITUDE, year:time, fl = LENGTH, 
         tagCode = TAG_CODE, broodYr = BROOD_YEAR, hatchery = HATCHERY,
         relSite = REL_SITE, reg = PROV_STATE, minRelJulian, maxRelJulian,
         pres = PRESERV)

write.csv(barkCWTTrim, here::here("data", "hargreaves", 
                                  "barkleyCWT_cleaned.csv"), row.names = FALSE)
