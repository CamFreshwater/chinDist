## Clean individual genetic data
# Data from west coast Vancouver Island commercial troll fishery
# Sep 14, 2019

library(tidyverse)
library(ggplot2)

datRaw <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbs_trim.csv"), 
                stringsAsFactors = FALSE)

# pull ID column, split and rename
id <- data.frame(t(matrix(unlist(strsplit(as.vector(datRaw$sampleID), 
                                          split = " ")),
                        nrow = 5, ncol = length(datRaw$sampleID)))) %>% 
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

dat <- cbind(id, datRaw) %>% 
  select(-sampleID, -abbYear) %>% 
  # for now remove ~150 fish that can't be assigned to an individual stat area
  filter(!statArea %in% c("Area123-124", "Area125-126", "Area126-127", 
                          "Area124_24")) %>% 
  select(-stockCode1, -stockCode2, -stockCode3, -stockCode4, -stockCode5) %>% 
  mutate(tempFish = paste(statArea, year, gear, jDay, fishNum, sep = "_"))


# collapse into regional groupings
stockKey <- read.csv(here::here("data", "southcoastStockKey.csv"), 
                     stringsAsFactors = FALSE) %>% 
  select(-Stock.Code, -Region1Code, -Region2Code, -Region3Code, 
         -FraserGroupCode)

gatherProbs <- dat %>% 
  select(tempFish, prob1, prob2, prob3, prob4, prob5) %>% 
  gather(., key = "probClass", value = "prob", -tempFish) %>% 
  arrange(tempFish)

gatherStocks <- dat %>% 
  select(tempFish, stock1, stock2, stock3, stock4, stock5) %>% 
  gather(., key = "stockClass", value = "Stock", -tempFish) %>% 
  left_join(., stockKey, by = "Stock") %>% 
  arrange(tempFish) %>% 
  cbind(., gatherProbs %>% select(-tempFish))

reg3 <- gatherStocks %>% 
  select(-Region2Name, -Region1Name) %>% 
  group_by(tempFish, Region3Name) %>% 
  dplyr::summarize(aggProb = sum(prob)) %>% 
  dplyr::arrange(tempFish, desc(aggProb)) %>% 
  right_join(dat %>% 
              select(tempFish, statArea, year, month, week, jDay, gear, fishNum,
                     date),
             .,
             by = "tempFish") %>% 
  mutate(Region = NA) %>% 
  dplyr::rename(regName = Region3Name)

reg3 <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                    "reg3RollUpCatchProb.csv")) %>% 
  select(-X)

## Code snippet to add regional identifiers
# fishSeq <- unique(reg3$tempFish)
# for(i in seq_along(fishSeq)) {
#   dum <- reg3 %>%
#     filter(tempFish == fishSeq[i])
#   for(j in 1:nrow(dum)) {
#     reg3[reg3$tempFish == fishSeq[i], ]$Region[j] <- paste("Region", j, 
#                                                            sep = " ")
#   }
# }
# write.csv(reg3, here::here("data", "gsiCatchData", "commTroll",
#                            "reg3RollUpCatchProb.csv"))


## Summarize samples and incorporate catch data
#monthly catches
# catch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
#                              "fosCatch.csv"),
#                   stringsAsFactors = FALSE) %>%
#   dplyr::group_by(MGMT_AREA, FISHING.YEAR, FISHING.MONTH) %>%
#   dplyr::summarize(catch = sum(CHINOOK_KEPT), boatDays = sum(VESSELS_OP)) %>%
#   dplyr::rename(statArea = "MGMT_AREA", year = "FISHING.YEAR",
#          month = "FISHING.MONTH") %>%
#   ungroup() %>%
#   mutate(statArea = as.character(statArea),
#          year = as.character(year))

#weekly catches
catch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                             "dailyCatch_WCVI.csv"),
                  stringsAsFactors = FALSE) %>% 
  dplyr::rename(statArea = area) %>% 
  mutate(date = as.Date(as.numeric(as.character(jDay)), 
                        origin = as.Date(paste(year, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
         week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d"))) %>% 
  select(-catchReg)

weeklyCatch <- catch %>%
  select(-jDay, - date, -sumCPUE) %>% 
  dplyr::group_by(statArea, year, month, week) %>% 
  dplyr::summarize(weeklyCatch = sum(catch),
            weeklyEffort = sum(boatDays))

reg3Catch <- reg3 %>% 
  mutate(statArea = as.numeric(as.character(statArea)),
         year = as.numeric(as.character(year))) %>% 
  full_join(., 
            weeklyCatch, 
            by = c("statArea", "year", "week", "month"))

## Distribution of samples and sampling effort
weekSamples <- dat %>% 
  group_by(statArea, year, month, week) %>% 
  tally() %>% 
  ungroup() %>%
  mutate(statArea = as.numeric(as.character(statArea)),
         year = as.numeric(as.character(year)),
         month = as.numeric(as.character(month)),
         week = as.numeric(as.character(week)))
# hist(weekSamples$n, breaks = 30)

summDat <- expand.grid(unique(dat$statArea), unique(dat$year), 
                       unique(dat$month), unique(dat$week)) %>%
  dplyr::rename("statArea" = Var1, "year" = Var2, "month" = Var3, 
                "week" = Var4) %>%
  mutate(statArea = as.numeric(as.character(statArea)),
         year = as.numeric(as.character(year))) %>% 
  full_join(., weekSamples, by = c("statArea", "year", "week", "month")) %>% 
  full_join(., weeklyCatch, by = c("statArea", "year", "week", "month")) %>% 
  replace_na(list(n = 0)) %>% 
  mutate(sampPpn = n / weeklyCatch) %>% 
  filter(statArea %in% unique(reg3$statArea)) #constrain to focal stat areas

ggplot(summDat, aes(x = week, y = sampPpn)) +
  geom_point() +
  facet_wrap(~statArea, scales = "free_y") +
  theme_bw()
ggplot(summDat, aes(x = week, y = n)) +
  geom_point() +
  facet_wrap(~statArea, scales = "free_y") +
  theme_bw()

## Sampling proportion exceeds 100% because at least some samples are from Taaq
# fishery; ideally include effort or at least catch from that sector
## Double check this...
dum2 <- summDat %>% 
  filter(sampPpn > 1)
write.csv(dum2, here::here("data", "gsiCatchData", "commTroll", 
                           "missingCatchData.csv"), row.names = FALSE)

summDat %>% 
  filter(statArea %in% c("24", "26"), 
         n > 0)
reg3 %>% 
  filter(statArea %in% c("24", "26"))
