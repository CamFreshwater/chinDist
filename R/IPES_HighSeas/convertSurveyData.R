## Reformat ship data
# Reformat data collected as .csv's during surveys to match datatable format
# of high seas salmon access database
# Oct. 9, 2019

library(tidyverse)


### Format to match BRIDGE Access table
bridgeDat <- read.csv(here::here("data", "highSeas", "importHighSeas", 
                    "bridge_shipData.csv"), stringsAsFactors = FALSE)
#convert header to all caps
# colnames(bridgeDat) <- stringr::str_to_upper(colnames(bridgeDat))

## Helper functions 
# Convert lat/longs to Access format
locConvert <- function(degreeIn, minuteIn) {
  paste(degreeIn, minuteIn, sep = " ") %>% 
    measurements::conv_unit(., from = "deg_dec_min", to = "dec_deg") %>% 
    as.numeric()
}

# Calculate distance
calcDist <- function(lat_a, lon_a, lat_b, lon_b) { 
  library(geosphere)
  distm(c(lon_a, lat_a), c(lon_b, lat_b), fun = distGeo) 
} 

# Calculate duration of set 
calcDuration <- function(startTime, stopTime) {
  start <- as.POSIXct(startTime,  format = "%H:%M")
  stop <- as.POSIXct(stopTime, format = "%H:%M")

  library(lubridate)
  interval(ymd_hms(start), ymd_hms(stop)) %>% 
    time_length(., "minute") / 60
}
  
# Generate output file
bridgeOut <- bridgeDat %>% 
  filter(Event_Type == "SET") %>% 
  mutate(STATION_ID = paste("BCSI", "2019125", Station, sep = "-"), 
         DATE_temp = as.POSIXct(paste(Year, Month, Day, sep = "-"), 
                           format = "%Y-%m-%d"),
         MONTH = DATE_temp %>% #upper case abbreviated
           lubridate::month(., label = TRUE, abbr = TRUE) %>% 
           stringr::str_to_upper(.),
         month_temp <- DATE_temp %>% #lower case unabbreviated
           lubridate::month(., label = TRUE),
         DATE = paste(Year, month_temp, Day, sep = "-"),
         JULIAN_DATE = lubridate::yday(DATE_temp),
         SHIP_TIME_ZONE = "PDT",
         START_LAT = locConvert(Start_Lat_Deg, Start_Lat_Min),
         START_LONG = locConvert(Start_Long_Deg, Start_Long_Min),
         END_LAT = locConvert(Stop_Lat_Deg, Stop_Lat_Min),
         END_LONG = locConvert(Stop_Long_Deg, Stop_Long_Min),
         DISTANCE = mapply(lat_a =  START_LAT, lon_a = START_LONG, #calc. dist.
                           lat_b = END_LAT, lon_b = END_LONG,
                           FUN = calcDist) / 1852, #division converts to naut. miles
         DUR = calcDuration(Start_Time, Stop_Time) %>% 
           round(., digits = 2),
         `SOG-KNTS` = Avg_Speed,
         HEADING = NA, #not sure how to calculate this and didn't have internet,
         START_BOT_DEPTH = "", #not recorded
         END_BOT_DEPTH = "", #not recorded
         `STAR/PORT` = "", #not recorded
         NET_OPENING_HEIGHT = "", #not recorded
         HEAD_DEPTH = Avg_Head_Depth,
         `ESTIMATED_CATCHES_BY_WT` = "") %>%  
  select(STATION_ID, EVENT = Event, STATION = Station, YEAR = Year, MONTH, 
         DAY = Day, DATE, JULIAN_DATE, START_TIME = Start_Time, SHIP_TIME_ZONE,
         START_LAT, START_LONG, END_LAT, END_LONG, DISTANCE, DUR, `SOG-KNTS`,
         START_BOT_DEPTH, END_BOT_DEPTH, `STAR/PORT`, 
         NET_OPENING_WIDTH = Doorspread, HEAD_DEPTH = Avg_Head_Depth, 
         WARP = Avg_Warp, `ESTIMATED_CATCHES_BY_WT`)


