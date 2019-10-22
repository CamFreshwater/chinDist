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
stockKey <- read.csv(here::here("data", "southcoastStockKey.csv"), 
                     stringsAsFactors = FALSE) %>% 
  select(stock = Stock, Region1Name, Region2Name, Region3Name)

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


#Merge id vector with original data frame and trim
dat <- cbind(id, datRaw) %>% 
  #calculate total summed probability for each sample
  group_by(szLineInfo) %>%
  mutate(totalProb = sum(dProb)) %>% 
  ungroup() %>% 
  mutate(adjProb = dProb / totalProb) %>% 
  select(-iRun, -iSample, -iYearMix, tempFish = szLineInfo, stock = szStock, 
         prob = dProb, adjProb, -totalProb, -szExclude, -iRegionId, 
         -szRegion) %>% 
  # for now remove ~150 fish that can't be assigned to an individual stat area
  filter(!statArea %in% c("Area123-124", "Area125-126", "Area126-127",
                          "Area124_24")) %>%
  left_join(., stockKey, by = "stock") #add southcoast regional groupings


#Roll up to regional aggregates (region 3 first)
reg3 <- dat %>% 
  group_by(tempFish, Region3Name) %>% 
  dplyr::summarize(aggProb = sum(adjProb)) %>% 
  dplyr::arrange(tempFish, desc(aggProb)) %>% 
  left_join(dat %>% 
               select(tempFish, statArea, year, month, week, jDay, gear, fishNum,
                      date),
             .,
             by = "tempFish") %>% 
  distinct() %>% 
  dplyr::rename(regName = Region3Name)

# write.csv(reg3, here::here("data", "gsiCatchData", "commTroll",
#                            "reg3RollUpCatchProb.csv"))

#Add weekly catches and sampling effort
weeklyCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                             "dailyCatch_WCVI.csv"),
                  stringsAsFactors = FALSE) %>% 
  dplyr::rename(statArea = area) %>% 
  mutate(statArea = as.character(statArea),
         date = as.Date(as.numeric(as.character(jDay)), 
                        origin = as.Date(paste(year, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
         week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d"))) %>% 
  select(-catchReg) %>%
  select(-jDay, - date, -sumCPUE) %>% 
  dplyr::group_by(statArea, year, month, week) %>% 
  dplyr::summarize(weeklyCatch = sum(catch),
                   weeklyEffort = sum(boatDays))
# write.csv(weeklyCatch, here::here("data", "gsiCatchData", "commTroll",
#                                    "weeklyCatch_WCVI.csv"))

weeklySamples <- dat %>% 
  group_by(statArea, year, month, week) %>% 
  summarize(nSampled = length(unique(fishNum))) %>% 
  # distinct() %>% 
  # tally(., name = "nSampled") %>% 
  ungroup() %>%
  mutate(statArea = as.character(statArea),
         year = as.numeric(as.character(year)),
         month = as.numeric(as.character(month)),
         week = as.numeric(as.character(week)))

reg3Catch <- reg3 %>% 
  mutate(statArea = as.character(statArea),
         year = as.numeric(as.character(year))) %>% 
  full_join(., 
            weeklyCatch, 
            by = c("statArea", "year", "week", "month")) %>% 
  left_join(., 
            weeklySamples,
            by = c("statArea", "year", "week", "month"))

# review sampling effort relative to catch
summDat <- weeklyCatch %>%
  full_join(.,
            weeklySamples,
            by = c("statArea", "year", "week", "month")) %>% 
  replace_na(list(nSampled = 0)) %>% 
  mutate(sampPpn = nSampled / weeklyCatch) %>% 
  # filter(statArea %in% unique(reg3$statArea)) %>% #constrain to focal stat areas
  distinct()

ggplot(summDat, aes(x = week, y = sampPpn)) +
  geom_point() +
  facet_wrap(~statArea, scales = "free_y") +
  theme_bw()
ggplot(summDat, aes(x = week, y = nSampled)) +
  geom_point() +
  facet_wrap(~statArea, scales = "free_y") +
  theme_bw()

## Sampling proportion exceeds 100% because at least some samples are from Taaq
# fishery; ideally include effort or at least catch from that sector
## Double check this...
dum2 <- summDat %>% 
  filter(sampPpn > 1) %>% 
  arrange(year)
write.csv(dum2, here::here("data", "gsiCatchData", "commTroll", 
                           "missingCatchData.csv"), row.names = FALSE)

