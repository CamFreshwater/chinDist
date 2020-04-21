## Clean individual genetic data - Ver 2
# Updated version of gsiIndProbsClean.R that uses new flat file with full
# probabilities (i.e. not trimmed to 5 stocks max) for WCVI commercial troll 
# fishery
# Oct 10, 2019

library(tidyverse)
library(ggplot2)

datRaw <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                              "wcviIndProbsLong_RAW.txt"), 
                   stringsAsFactors = FALSE)

# Big chunk of code to separate ID variable into meaningful individual vectors
id <- datRaw$szLineInfo %>% 
  as.vector() %>% 
  strsplit(., split = " ") %>% 
  unlist() %>% 
  matrix(., nrow = 5, ncol = length(datRaw$szLineInfo)) %>%
  t() %>% 
  data.frame() %>% 
  rename("statArea" = X1, "year" = X2, "gear" = X3, "jDay" = X4,
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

# Export example GSI data to Wilf 
# dat[1:1000, ] %>% 
#   select(flatFileID:prob, stock, region = Region1Name, statArea, gear, 
#          fish_number = fishNum, 
#          date, year, month, yday = jDay) %>% 
#   write.csv(., here::here("data", "gsiCatchData", "commTroll",
#                           "example_gsi_input.csv"), row.names = FALSE)

## Export list of stocks to be passed to makeFullStockKey script in
# stockKey repo 
# stks_out <- dat %>%
#   select(stock, Region1Name) %>%
#   distinct()
# saveRDS(stks_out, here::here("data", "stockKeys", "wcviTrollStocks.rds"))

# stock key generated in stockKey repo
stockKey <- readRDS(here::here("data", "stockKeys", "finalStockList_Mar2020.rds"))

dat2 <- dat %>%
  select(-Region1Name) %>%
  left_join(., stockKey, by = c("stock")) %>% 
  arrange(statArea, year, month, jDay, fishNum) %>% 
  rename(id = flatFileID)
# saveRDS(dat2, here::here("data", "gsiCatchData", "commTroll",
#                          "wcviIndProbsLong.rds"))
dat2 <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                           "wcviIndProbsLong.rds"))


## Roll up to regional aggregates ----------------------------------------------

#cleaning function to filter out non-dominant assignments below threshold
clean_dat <- function(dat, threshold = 0.75) {
  dat %>% 
    group_by(flatFileID) %>% 
    mutate(max_assignment = max(aggProb)) %>% 
    # Remove samples where top stock ID is less than 75% probability
    filter(!aggProb < max_assignment, 
           !max_assignment < threshold) %>% 
    ungroup() %>% 
    distinct() %>% 
    mutate(
      season_c = case_when(
        month %in% c("12", "1", "2") ~ "w",
        month %in% c("3", "4", "5") ~ "sp",
        month %in% c("6", "7", "8") ~ "su",
        month %in% c("9", "10", "11") ~ "f"
      ),
      season = fct_relevel(season_c, "sp", "su", "f", "w"),
      month_n = as.numeric(month),
      month = as.factor(month_n),
      year =  as.factor(year),
      pres = 1,
      area_n = as.numeric(as.character(statArea)),
      catchReg = case_when(
        area_n < 125 & area_n > 27 ~ "SWVI",
        area_n < 25 ~ "SWVI",
        TRUE ~ "NWVI"
      ),
      catchReg = as.factor(catchReg), 
      statArea = as.factor(statArea)) 
}

# Region 3 first (i.e. large regional aggregats)
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
  dplyr::rename(regName = Region3Name) %>% 
  # aggregate based on likely stocks of interest
  mutate(regName = fct_recode(regName, ECVI = "SOG"),
         regName = as.character(regName),
         regName = case_when(
           regName %in% c("Columbia", "Snake") ~ "Columbia",
           regName %in% c("Coastal Washington", "Washington Coast",
                          "Alaska South SE", "North/Central BC", "SOG", 
                          "Oregon/California", "ECVI", "WCVI") ~ "Other",
           TRUE ~ regName
         ),
         regName = as.factor(abbreviate(regName, minlength = 5))
  ) %>% 
  #remove non-dom assignments based on above
  clean_dat(threshold = 0.75)

saveRDS(reg3, here::here("data", "gsiCatchData", "commTroll",
                           "reg3RollUpCatchProb.RDS"))

#key for adding aggregates to below 
reg3_key <- reg3 %>%
  select(flatFileID, aggName = regName)

# Region 1 next (approximately equivalent to PSC groupings eg MUFR)
reg1 <- dat2 %>% 
  group_by(flatFileID, Region1Name) %>% 
  dplyr::summarize(aggProb = sum(adjProb)) %>% 
  dplyr::arrange(flatFileID, desc(aggProb)) %>% 
  left_join(dat %>% 
              select(flatFileID, statArea, year, month, week, jDay, gear, 
                     fishNum, date),
            .,
            by = "flatFileID") %>% 
  distinct() %>% 
  dplyr::rename(pscName = Region1Name) %>% 
  clean_dat() %>% 
  left_join(.,
            reg3_key,
            by = "flatFileID")

saveRDS(reg1, here::here("data", "gsiCatchData", "commTroll",
                         "reg1RollUpCatchProb.RDS"))

tt <- fr_reg1B %>% 
  filter(aggName == "FrsrR")
table(tt$pscName, tt$month_n)

# Modified region 1 with Fraser focus
fr_reg1 <- reg1 %>%
  mutate(
    regName = case_when(
      aggName == "FrsrR" & grepl("TH", pscName) ~ "Thomp-Early",
      pscName %in% c("LWFR-Su", "LWFR-Sp", "MUFR", "UPFR") ~ "FR-Early",
      pscName == "LWFR-F" ~ "LWFR-Late",
      aggName == "FrsrR" ~ pscName,
      TRUE ~ "Other"
    )
  ) 
saveRDS(fr_reg1, here::here("data", "gsiCatchData", "commTroll",
                            "reg1RollUpCatchProb_Fraser.RDS"))
fr_reg1 <- saveRDS(here::here("data", "gsiCatchData", "commTroll",
                            "reg1RollUpCatchProb_Fraser.RDS"))
# alternative option w/ finer resolution 
fr_reg1B <- reg1 %>%
  mutate(
    regName = case_when(
      pscName == "SOTH" ~ "South Thomp.",
      pscName %in% c("LWTH", "NOTH") ~ "Thomp. 1.x",
      pscName %in% c("LWFR-Su", "LWFR-Sp", "MUFR", "UPFR") ~ "FR 1.x",
      pscName == "LWFR-F" ~ "FR-Fall",
      aggName == "FrsrR" ~ pscName,
      TRUE ~ "Other"
    )
  ) 
saveRDS(fr_reg1B, here::here("data", "gsiCatchData", "commTroll",
                            "reg1RollUpCatchProb_FraserB.RDS"))
table(fr_reg1$regName, fr_reg1$month_n)
table(fr_reg1B$regName, fr_reg1B$month_n)


# Region 2 next (intermediate to 1 and 3)
reg2 <- dat2 %>% 
  group_by(flatFileID, Region2Name) %>% 
  dplyr::summarize(aggProb = sum(adjProb)) %>% 
  dplyr::arrange(flatFileID, desc(aggProb)) %>% 
  left_join(dat %>% 
              select(flatFileID, statArea, year, month, week, jDay, gear, 
                     fishNum, date),
            .,
            by = "flatFileID") %>% 
  distinct() %>% 
  dplyr::rename(pscName = Region2Name) %>% 
  clean_dat() %>% 
  left_join(.,
            reg3_key,
            by = "flatFileID")

saveRDS(reg2, here::here("data", "gsiCatchData", "commTroll",
                         "reg2RollUpCatchProb.RDS"))


# Compare to catch data --------------------------------------------------------
#Add weekly catches and sampling effort
dailyCatch <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                    "dailyCatch_WCVI.rds")) %>% 
  dplyr::rename(statArea = area) %>% 
  mutate(statArea = as.character(statArea),
         date = as.Date(as.numeric(as.character(jDay)), 
                        origin = as.Date(paste(year, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
         week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d")))
weeklyCatch <- dailyCatch %>% 
  select(-catchReg) %>%
  select(-jDay, - date, -cpue) %>% 
  dplyr::group_by(statArea, year, month, week) %>% 
  dplyr::summarize(weeklyCatch = sum(catch),
                   weeklyEffort = sum(boatDays))
write.csv(weeklyCatch, here::here("data", "gsiCatchData", "commTroll",
                                   "weeklyCatch_WCVI.csv"), row.names = FALSE)


# review sampling effort relative to catch
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