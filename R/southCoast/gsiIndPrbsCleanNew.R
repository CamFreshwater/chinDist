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
# stock key generated in juvenile salmon index repo
stockKey <- readRDS(here::here("data", "stockKeys", "finalStockList.rds")) 

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
  mutate(jDay = as.numeric(as.character(jDay)),
         statArea = 
           case_when(
              statArea %in% c("Area023", "Area23", "Area_23") ~ "23",
              statArea %in% c("Area123", "Area123SWVI", "Area123Comm") ~ "123",
              statArea %in% c("Area124", "Area124SWVI", "Area124Comm") ~ "124",
              statArea %in% c("Area125", "Area125NWVI") ~ "125",
              statArea %in% c("Area126", "Area126NWVI") ~ "126",
              statArea %in% c("Area127", "Area127NWVI") ~ "127",
              statArea %in% c("Area026") ~ "26",
              statArea %in% c("Area24", "Area_24", "Area_24xgill") ~ "24",
              TRUE ~ as.character(statArea)
              ),
         abbYear = sapply(strsplit(as.character(year), '[()]'), 
                     function(x) (x)[2]),
         year = paste("20", abbYear, sep = ""),
         #adjust sampling day to correct for errors by genetics lab and 
         jDay = 
           case_when(
             statArea == "126" & year == "2012" & jDay > 48 &
               jDay < 116 ~ 48,
             statArea == "126" & year == "2012" & jDay > 116 &
               jDay < 121 ~ 116,
             TRUE ~ jDay),
         date = as.Date(as.numeric(as.character(jDay)),
                        origin = as.Date(paste(year, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
         week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d")),
         fishNum = as.numeric(as.character(fishNum))) %>% 
  # adjust sampling week
  # for specific strata based on when samples were landed relative to julian 
  # date of when fishing occurred (based on reviewing FOS database and pers. 
  # comm. Lee Kearey - South Coast)
  mutate(adjWeek = case_when(
    statArea == "23" & month == "3" & year == "2013" ~ week - 1,
    statArea == "24" & month == "5" & year == "2013" ~ week - 1,
    statArea == "26" & month == "2" & year == "2" ~ week - 1,
    statArea == "123" & month == "6" & year == "2007" ~ week - 1,
    statArea == "123" & month == "6" & year == "2008" ~ week - 1,
    statArea == "125" & month == "6" & year == "2007" ~ week - 1,
    statArea == "125" & month == "10" & year == "2007" ~ week - 1,
    statArea == "125" & month == "5" & year == "2014"~ week - 1,
    statArea == "126" & month == "3" & year == "2007"~ week - 1,
    statArea == "126" & month == "2" & year == "2012"~ week - 1,
    statArea == "126" & month == "3" & year == "2013"~ week - 1,
    statArea == "126" & month == "4" & year == "2013"~ week - 1,
    statArea == "126" & month == "3" & year == "2015"~ week - 1,
    statArea == "127" & month == "6" & year == "2007"~ week - 1,
    statArea == "127" & month == "4" & year == "2014"~ week - 1,
    statArea == "23" & month == "6" & year == "2007"~ week - 1,
    TRUE ~ week)
  ) %>% 
  #shuffle adjusted
  rename(week = adjWeek, unadjWeek = week)

#Merge id vector with original data frame and trim
dat <- cbind(id, datRaw) %>% 
  #calculate total summed probability for each sample
  group_by(szLineInfo) %>%
  mutate(stock = toupper(szStock),
         totalProb = sum(dProb)) %>% 
  ungroup() %>% 
  mutate(adjProb = dProb / totalProb) %>% 
  select(-iRun, -iSample, -iYearMix, flatFileID = szLineInfo, stock, -szStock,
         prob = dProb, adjProb, -totalProb, -szExclude, -iRegionId, -unadjWeek,
         Region1Name = szRegion) %>% 
  # for now remove ~150 fish that can't be assigned to an individual stat area;
  # eventually could assign based on where majority of effort occurred
  filter(!statArea %in% c("Area123-124", "Area125-126", "Area126-127",
                          "Area124_24")) 

## Export list of stocks to be passed to makeFullStockKey script in
# juvenile-salmon-index repo 
# stks_out <- dat %>%
#   select(stock, Region1Name) %>%
#   distinct()
# saveRDS(stks_out, here::here("data", "stockKeys", "wcviTrollStocks.rds"))

dat2 <- dat %>% 
  select(-Region1Name) %>% 
  left_join(., stockKey, by = c("stock"))
saveRDS(dat2, here::here("data", "gsiCatchData", "commTroll",
                         "wcviIndProbsLong_CLEAN.rds"))

dailySamples <- dat2 %>% 
  group_by(statArea, year, month, jDay) %>% 
  mutate(nSampled = length(unique(fishNum))) %>% 
  ungroup() %>%
  select(statArea, year, month, week, jDay, nSampled) %>% 
  distinct() %>% 
  mutate(year = as.numeric(year))
weeklySamples <- dailySamples %>% 
  group_by(statArea, year, month, week) %>% 
  summarize(nSampled = sum(nSampled)) %>% 
  ungroup()

#Roll up to regional aggregates (region 3 first)
reg3 <- dat2 %>% 
  group_by(flatFileID, Region3Name) %>% 
  dplyr::summarize(aggProb = sum(adjProb)) %>% 
  dplyr::arrange(flatFileID, desc(aggProb)) %>% 
  left_join(dat %>% 
               select(flatFileID, statArea, year, month, week, jDay, gear, 
                      fishNum, date),
             .,
             by = "flatFileID") %>% 
  distinct() %>% 
  dplyr::rename(regName = Region3Name)

write.csv(reg3, here::here("data", "gsiCatchData", "commTroll",
                           "reg3RollUpCatchProb.csv"), row.names = FALSE)

#Add weekly catches and sampling effort
dailyCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                    "dailyCatch_WCVI.csv"),
         stringsAsFactors = FALSE) %>% 
  dplyr::rename(statArea = area) %>% 
  mutate(statArea = as.character(statArea),
         date = as.Date(as.numeric(as.character(jDay)), 
                        origin = as.Date(paste(year, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
         week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d")))
weeklyCatch <- dailyCatch %>% 
  select(-catchReg) %>%
  select(-jDay, - date, -sumCPUE) %>% 
  dplyr::group_by(statArea, year, month, week) %>% 
  dplyr::summarize(weeklyCatch = sum(catch),
                   weeklyEffort = sum(boatDays))
write.csv(weeklyCatch, here::here("data", "gsiCatchData", "commTroll",
                                   "weeklyCatch_WCVI.csv"), row.names = FALSE)


# review sampling effort relative to catch
summDat <- weeklyCatch %>%
  full_join(.,
            # adjWeeklySamps,
            weeklySamples,
            by = c("statArea", "year", "week", "month")) %>% 
  replace_na(list(nSampled = 0, weeklyCatch = 0, weeklyEffort = 0)) %>% 
  mutate(sampPpn = nSampled / weeklyCatch) %>% 
  # filter(statArea %in% unique(reg3$statArea)) %>% #constrain to focal stat areas
  distinct()

## Sampling proportion exceeds 100% because at least some samples are from Taaq
# fishery; ideally include effort or at least catch from that sector
## Double check this...
dum2 <- summDat %>% 
  filter(sampPpn > 1) %>% 
  arrange(statArea, year) 

dailyCatch %>% 
  filter(statArea == "24", month == "6", year == "2013") %>% 
  select(statArea:jDay, week, catch, date)
dailySamples %>% 
  filter(statArea == "24", month == "6", year == "2013")

# weeklyCatch %>% 
#   filter(statArea == "125", month == "9", year == "2012") %>% 
#   select(-weeklyEffort)
# weeklySamples %>% 
#   filter(statArea == "125", month == "9", year == "2012")

write.csv(dum2, here::here("data", "gsiCatchData", "commTroll", 
                           "missingCatchData.csv"), row.names = FALSE)

ggplot(summDat, aes(x = week, y = sampPpn)) +
  geom_point() +
  facet_wrap(~statArea, scales = "free_y") +
  theme_bw()
ggplot(summDat, aes(x = week, y = nSampled)) +
  geom_point() +
  facet_wrap(~statArea, scales = "free_y") +
  theme_bw()