## Clean individual genetic data - Ver 2
# Updated version of gsiIndProbsClean.R that uses new flat file with full
# probabilities (i.e. not trimmed to 5 stocks max) for WCVI commercial troll 
# fishery
# Oct 10, 2019

library(tidyverse)
library(ggplot2)

datRaw <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                              "wcviIndProbsLong.txt"), 
                   stringsAsFactors = FALSE)

# Big chunk of code to separate ID variable into meaningful individual vectors
id <- datRaw$szLineInfo %>% 
  as.vector() %>% 
  strsplit(., split = " ") %>% 
  unlist() %>% 
  matrix(., nrow = 5, ncol = length(datRaw$szLineInfo)) %>%
  t() %>% 
  data.frame() %>% 
  dplyr::rename("statArea" = X1, "year" = X2, "gear" = X3, "jDay" = X4, 
                "fishNum" = X5) %>% 
  mutate(statArea = case_when(
    statArea %in% c("Area023", "Area23", "Area_23") ~ "23",
    statArea %in% c("Area123", "Area123SWVI", "Area123Comm") ~ "123",
    statArea %in% c("Area124", "Area124SWVI", "Area124Comm") ~ "124",
    statArea %in% c("Area125", "Area125NWVI") ~ "125",
    statArea %in% c("Area126", "Area126NWVI") ~ "126",
    statArea %in% c("Area127", "Area127NWVI") ~ "127",
    statArea %in% c("Area026") ~ "26",
    statArea %in% c("Area24", "Area_24", "Area_24xgill") ~ "24",
    TRUE ~ as.character(statArea)),
    abbYear = sapply(strsplit(as.character(year), '[()]'), 
                     function(x) (x)[2]),
    year = paste("20", abbYear, sep = ""),
    date = as.Date(as.numeric(as.character(jDay)), 
                   origin = as.Date(paste(year, "01", "01", sep = "-"))),
    month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
    week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d")))

# Merge id vector with original data frame and trim
dat <- cbind(id, datRaw) %>% 
  select(-iRun, -iSample, -iYearMix, -szLineInfo, -szStock, prob = dProb, 
         -szExclude, -iRegionId, genRegion = szRegion) 

%>% 
  # for now remove ~150 fish that can't be assigned to an individual stat area
  filter(!statArea %in% c("Area123-124", "Area125-126", "Area126-127", 
                          "Area124_24")) %>% 
  select(-stockCode1, -stockCode2, -stockCode3, -stockCode4, -stockCode5) %>% 
  mutate(tempFish = paste(statArea, year, gear, jDay, fishNum, sep = "_"))
