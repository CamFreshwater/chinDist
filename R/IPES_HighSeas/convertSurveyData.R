## Reformat ship data
# Reformat data collected as .csv's during surveys to match datatable format
# of high seas salmon access database
# Oct. 9, 2019

library(tidyverse)

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

# Filter count 
calcCount <- function(speciesCode)
  haulDat %>% 
  filter(speciesCode)


## Data recorded on board
bridgeDat <- read.csv(here::here("data", "highSeas", "importHighSeas", 
                                 "bridge_shipData.csv"), 
                      stringsAsFactors = FALSE)
indDat <- read.csv(here::here("data", "highSeas", "importHighSeas", 
                              "ind_shipData.csv"), stringsAsFactors = FALSE)
haulDat <- read.csv(here::here("data", "highSeas", "importHighSeas",
                               "haulCard_shipData.csv"), 
                    stringsAsFactors = FALSE) %>% 
  filter(CRUISE.NUMBER == "2019-125")



############
##### Format to match BRIDGE Access table

# Calculate summary catch for coho and chinook size classes
# dummy dataset with 0 catches for each size bin to ensure all are spread
dum <- data.frame(
  STATION = "DUM", 
  FISH.NUMBER = rep(seq(from = 1, to = 11, by = 1), times = 2),
  SPECIES = c(rep("COHO", times = 11), rep("CHINOOK", times = 11)),
  FL = rep(seq(from = 98, to = 1098, by = 100), times = 2)
)

piscDat <- indDat %>% 
  filter(SPECIES %in% c("COHO", "CHINOOK")) %>% 
  select(STATION, FISH.NUMBER, SPECIES, FL = "FL..mm.") %>% 
  rbind(., dum) %>% #add dummy dataset
  mutate(flCat = 
           case_when(
             FL < 100 ~ "<100MM",
             100 < FL & FL < 199 ~ "100-199MM",
             199 < FL & FL < 299 ~ "200-299MM",
             299 < FL & FL < 399 ~ "300-399MM",
             399 < FL & FL < 499 ~ "400-499MM",
             499 < FL & FL < 599 ~ "500-599MM",
             599 < FL & FL < 699 ~ "600-699MM",
             699 < FL & FL < 799 ~ "700-799MM",
             799 < FL & FL < 899 ~ "800-899MM",
             899 < FL & FL < 999 ~ "900-999MM",
             999 < FL & FL < 1099 ~ "1000-1099MM"),
         abbSp = 
           case_when(
             SPECIES == "COHO" ~ "CO",
             SPECIES == "CHINOOK" ~ "CK"),
         spCat = paste(abbSp, flCat, sep = "_")) %>% 
  group_by(STATION, spCat) %>% 
  summarize(count = length(unique(FISH.NUMBER))) %>% 
  spread(., spCat, count) %>% 
  filter(!STATION == "DUM") %>% 
  select(-`CO_1000-1099MM`)

# Calculate summary catch for all salmon species
summCatch <- haulDat %>%
  mutate(sp_t = case_when(
    SPECIES.CODE == "108J" ~ "PK_JUV",
    SPECIES.CODE == "124J" ~ "CK_JUV",
    SPECIES.CODE == "124A" ~ "CK_ADULT",
    SPECIES.CODE == "112J" ~ "CM_JUV",
    SPECIES.CODE == "112A" ~ "CM_ADULT",
    SPECIES.CODE == "115J" ~ "CO_JUV",
    SPECIES.CODE == "118J" ~ "SE_JUV",
    TRUE ~ "X"
  )) %>% 
  filter(!sp_t == "X") %>% 
  select(STATION, COUNT, sp_t) %>% 
  spread(., sp_t, COUNT) %>% 
  mutate(CO_ADULT = NA, PK_ADULT = NA, SE_ADULT = NA, SD_JUV = NA, 
         SD_ADULT = NA,
         CK = sum(CK_JUV + CK_ADULT)) %>% 
  left_join(., piscDat, by = "STATION")

# Quick loop to replace NAs with 0s (replace_na is cumbersome with so many columns)
for(i in 1:nrow(summCatch)) {
  seqNA <- which(is.na(summCatch[i, ]))
  summCatch[i, seqNA] <- 0
}

# Generate output file
bridgeTemp <- bridgeDat %>% 
  filter(Event_Type == "SET") %>% 
  mutate(STATION_ID = paste("BCSI", "2019125", Station, sep = "-"), 
         STATION = Station,
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
         `ESTIMATED_CATCHES_BY_WT` = "",
         GEAR_TYPE = "", #not recorded,
         COMMENTS = Comments)

bridgeOut <- bridgeTemp %>%  
  left_join(., summCatch, by = "STATION") %>% 
  select(STATION_ID, EVENT = Event, STATION, YEAR = Year, MONTH,
         DAY = Day, DATE, JULIAN_DATE, START_TIME = Start_Time, SHIP_TIME_ZONE,
         START_LAT, START_LONG, END_LAT, END_LONG, DISTANCE, DUR, `SOG-KNTS`,
         START_BOT_DEPTH, END_BOT_DEPTH, `STAR/PORT`,
         NET_OPENING_WIDTH = Doorspread, HEAD_DEPTH = Avg_Head_Depth,
         WARP = Avg_Warp, `ESTIMATED_CATCHES_BY_WT`, PK_JUV, PK_ADULT, CM_JUV,
         CM_ADULT, SE_JUV, SE_ADULT, CO_JUV, CO_ADULT, `CO_<100MM`,
         `CO_100-199MM`, `CO_200-299MM`, `CO_300-399MM`,`CO_400-499MM`,
         `CO_500-599MM`, `CO_600-699MM`, `CO_700-799MM`,`CO_800-899MM`,
         `CO_900-999MM`, CK, CK_JUV, CK_ADULT, `CK_100-199MM`, `CK_1000-1099MM`,
         `CK_200-299MM`, `CK_300-399MM`, `CK_400-499MM`, `CK_500-599MM`, 
         `CK_600-699MM`, `CK_700-799MM`, `CK_800-899MM`, `CK_900-999MM`,
         SD_JUV, SD_ADULT, GEAR_TYPE, COMMENTS)

write.csv(bridgeOut, here::here("data", "highSeas", "importHighSeas",
                                "bridgeToImport.csv"), row.names = FALSE)


############
##### Format to match BIOLOGICAL_JUNCTION Access table
juncOutTemp <- indDat %>% 
  mutate(STATION_ID = paste("BCSI", "2019125", STATION, sep = "-"),
         spCode = case_when(
           SPECIES == "PINK" ~ "108",
           SPECIES == "CHUM" ~ "112",
           SPECIES == "COHO" ~ "115",
           SPECIES == "SOCKEYE" ~ "118",
           SPECIES == "CHINOOK" ~ "124",
           TRUE ~ "NA"),
         CRUISE_STATION_ID_SPECIES = paste(STATION_ID, spCode, sep = "-"),
         FISH_NUMBER = sprintf("%03d", FISH.NUMBER) %>%
           paste(CRUISE_STATION_ID_SPECIES, ., sep = "-"),
         TM_CODE = "" #thermal mark code (not relevant)
         )

juncOut <- juncOutTemp %>% 
  filter(!spCode == "NA") %>% 
  select(STATION_ID, CRUISE_STATION_ID_SPECIES, FISH_NUMBER, CWT, TM_CODE)
  
write.csv(juncOut, here::here("data", "highSeas", "importHighSeas",
                              "biologicalJunctionToImport.csv"), 
          row.names = FALSE)


############
##### Format to match BIOLOGICAL Access table
samplingDate <- bridgeTemp %>% 
  dplyr::select(STATION_ID, DATE_temp, Day, Year) %>% 
  mutate(wDay = wday(DATE_temp, label = TRUE, abbr = FALSE),
         monthDay = DATE_temp %>% #upper case abbreviated
           lubridate::month(., label = TRUE, abbr = TRUE) %>% 
           paste(., Day, sep = " "),
         DATE_CAPTURED = paste(wDay, monthDay, Year, sep = ", ")) %>% 
  dplyr::select(STATION_ID, DATE_CAPTURED)

biolOut <- juncOutTemp %>% 
  left_join(., samplingDate, by = "STATION_ID") %>% 
  mutate(SHIP_SL = "",
         SHIP_TL = "",
         SHIP_ML = "",
         CWT_AGE = "",
         GONAD_WT = "",
         LICE_PREVALENCE = "",
         TEMPORARY_TAG = "") %>% 
  dplyr::select(FISH_NUMBER, DATE_CAPTURED, SPECIES, SHIP_FL = FL..mm.,
                SHIP_SL, SHIP_TL, SHIP_ML, SHIP_WT = WT..mm., 
                OTOLITH = OTOLITH.., CWT_AGE, GONAD_WT, FIN_CL = FIN.CLIP, SEX,
                SAMPLE, COMMENTS, TEMPORARY_TAG)
  
write.csv(biolOut, here::here("data", "highSeas", "importHighSeas",
                              "biologicalToImport.csv"), 
          row.names = FALSE)
  