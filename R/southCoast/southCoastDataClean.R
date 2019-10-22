## Clean data from southcoast GSI Access database
# Provided by Bryan Rusch Aug 19, 2019
# Note: the tables within the access database are not joined, limiting the
# utility of multitable queries; to ensure repeatability the majority of data
# cleaning occurs after tables have been imported to R

library(RODBC); library(tidyverse); library(ggplot2)

# Access database saved locally to work computer
file_path <- "C:/Users/FRESHWATERC/Documents/chinook/southCoastDatabase/southCoastGSICatch.accdb"
con <- odbcConnectAccess2007(file_path)

# Check out table names
sqlTables(con, tableType = "TABLE")$TABLE_NAME

## Import various datatables
# all gsi samples
sampInvQry <- "SELECT *
               FROM [Sample Inventory]"
sampInv <- sqlQuery(con, sampInvQry) %>% 
  select(-F14, -F15, -F16) %>% 
  rename_all(list(~make.names(.))) #get rid of spaces in col names
trimSampInv <- sampInv %>% 
  select(Sample.ID, Vial.ID, Year, Fishery, Catch.Region, Catch.Area, 
         Month.Code, Month.Name, Lab.Reported.N)

# stock composition data
stockCompQry <- "SELECT *
                 FROM [Stock Level Results] 
                 INNER JOIN [Stocks with Region Codes]
                 ON [Stock Level Results].[Stock Code] = 
                  [Stocks with Region Codes].[Stock Code];"
stockComp <- sqlQuery(con, stockCompQry) %>% 
  rename_all(list(~make.names(.))) %>% 
  left_join(., trimSampInv, by = "Sample.ID") %>% 
  select(Sample.ID:SD, Region1Code:Lab.Reported.N) %>% 
  rename(catchReg = Catch.Region, month = Month.Code, year = Year) %>% 
  mutate(catchReg = as.character(catchReg)) %>% 
  filter(Fishery == "Area G Com") 
yrs <- unique(stockComp$year)

# confirm one sample collected per catch region, year, and month
# dum <- stockComp %>%
#   select(-c(Stock.Code:Vial.ID)) %>%
#   distinct() %>%
#   group_by(Catch.Region, Year, Month.Name) %>%
#   tally()

# region key to reference stocks at different roll ups
reg1Qry <- "SELECT *
            FROM [Region 1 Stock Names];"
reg1 <- sqlQuery(con, reg1Qry) %>% 
  select(-ID)
reg3Qry <- "SELECT *
            FROM [Region 3 Stock Names];"
reg3 <- sqlQuery(con, reg3Qry) %>% 
  select(-ID)
# frNameQry <- "SELECT *
#               FROM [Fraser Stock Grouping Names]"
reg2Qry <- "SELECT *
           FROM [Stocks with Region Codes]
           LEFT JOIN [Region 2 Stock Names]
           ON [Stocks with Region Codes].Region2Code = 
              [Region 2 Stock Names].Region2Code;"
regKey <- sqlQuery(con, reg2Qry) %>%
  rename_all(list(~make.names(.))) %>%
  left_join(., reg1, by = "Region1Code") %>% 
  left_join(., reg3, by = "Region3Code") %>% 
  mutate(Stock = Stock.Name) %>% 
  select(Stock.Code, Stock, Region2Name, Region1Name, Region3Name, 
         Region1Code:FraserGroupCode)
head(regKey)

trimRegKey <- regKey %>% 
  mutate(Stock = as.character(Stock)) %>% 
  select(Stock, Region1Name, Region2Name, Region3Name)

# write.csv(regKey, here::here("data", "southcoastStockKey.csv"), row.names = F)

# all catch/effort data
catchQry <- "SELECT [Area G FOS Catch Estimates].ESTIMATE_TYPE, 
                    [Area G FOS Catch Estimates].LICENCE_AREA, 
                    [Area G FOS Catch Estimates].OPNG_CAT, 
                    [Area G FOS Catch Estimates].OPNG_DESC, 
                    [Area G FOS Catch Estimates].TARGETS_CHINOOK, 
                    [Area G FOS Catch Estimates].STAT_WEEK, 
                    [Area G FOS Catch Estimates].FISHING_DATE, 
                    Month([FISHING_DATE]) AS [FISHING MONTH], 
                    Year([FISHING_DATE]) AS [FISHING YEAR], 
                    [Area G FOS Catch Estimates].MGMT_AREA, 
                    [Area G FOS Catch Estimates].AREA_NAME, 
                    [Management Regions].CATCH_REGION, 
                    [Area G FOS Catch Estimates].HRS_OPEN, 
                    [Area G FOS Catch Estimates].VESSELS_OP, 
                    [Area G FOS Catch Estimates].CHINOOK_KEPT, 
                    [Area G FOS Catch Estimates].CHINOOK_RELD, 
                    [Area G FOS Catch Estimates].COMMENTS
               FROM [Area G FOS Catch Estimates] 
               INNER JOIN [Management Regions] 
               ON [Area G FOS Catch Estimates].MGMT_AREA = 
                  [Management Regions].[Managment Area]
               WHERE ((([Area G FOS Catch Estimates].TARGETS_CHINOOK)='Yes'));
"
areaCatch <- sqlQuery(con, catchQry) %>% 
  rename_all(list(~make.names(.))) %>%
  filter(FISHING.YEAR %in% yrs,
         !VESSELS_OP > 200) #remove one entry with unrealistically high effort
# write.csv(areaCatch, here::here("data", "gsiCatchData", "commTroll", 
#                                 "fosCatch.csv"))

areaCatch %>% 
  filter(FISHING.YEAR == "2011",
         FISHING.MONTH == "9",
         AREA_NAME == "management area 123") 

areaCatch %>% 
  filter(grepl("CN DNA Sampling", OPNG_DESC))

#summary of how catch/effort data is distributed through time
# fosSumm <- areaCatch %>%
#   group_by(CATCH_REGION, FISHING.MONTH, FISHING.YEAR, MGMT_AREA) %>%
#   tally(name = "daysWithData")
# 
# write.csv(fosSumm, here::here("data", "gsiCatchData", "commTroll", 
#                               "fosSummary.csv"))

#summarize catches by area, month and year to match south coast GSI data
sumCatch <- areaCatch %>% 
  group_by(CATCH_REGION, FISHING.MONTH, FISHING.YEAR) %>% 
  mutate(sumCatch = 
         #normally only KEPT chinook contribute to GSI samples, but certain
         #fisheries include live sampling where CHINOOK_RELD = catch
          case_when(
            grepl("CN DNA Sampling", OPNG_DESC) ~ sum(CHINOOK_RELD),
            TRUE ~ sum(CHINOOK_KEPT)),
         sumEffort = sum(VESSELS_OP)) %>% 
  select(catchReg = CATCH_REGION, month = FISHING.MONTH, year = FISHING.YEAR, 
         sumCatch, sumEffort) %>% 
  ungroup() %>% 
  mutate(catchReg = as.character(catchReg))

# # catch/effort coverage
# ggplot(sumCatch, aes(x = month, y = sumCatch, color = catchReg)) +
#   geom_line() +
#   facet_wrap(~year)
# ggplot(sumCatch, aes(x = month, y = sumEffort, color = catchReg)) +
#   geom_line() +
#   facet_wrap(~year)

# combine monthly summed catches w/ monthly GSI to calculate stock specific 
# catches (divide by 100 since they're percentages currently)
stockCatch <- stockComp %>% 
  left_join(., sumCatch, by = c("catchReg", "month", "year")) %>% 
  mutate(estCatch = (Estimate / 100) * sumCatch,
         varCatch = ((SD / 100) * sumCatch)^2,
         Stock = as.character(Stock),
         samplePpn = Lab.Reported.N / sumCatch) %>% 
  left_join(., trimRegKey, by = "Stock") %>% 
  select(catchReg, month, year, sumEffort, labN = Lab.Reported.N, samplePpn,
         estCatch, varCatch, Stock, Region1Name, Region2Name, Region3Name,
         Region1Code:FraserGroupCode)

dd <- stockCatch %>% 
  filter(is.na(samplePpn))

ee <- areaCatch %>% 
  filter(FISHING.YEAR == "2013", 
  FISHING.MONTH %in% c("6","8"))

write.csv(stockCatch, here::here("data", "gsiCatchData", "commTroll", 
                                 "stockCatch_WCVI.csv"), row.names = FALSE)


## Generate aggregate catch by Julian Day to match individual data from genetics
# lab
areaCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                                 "fosCatch.csv"))

dailyCatch <- areaCatch %>% 
  mutate(jDay = lubridate::yday(as.POSIXlt(FISHING_DATE, 
                                           format="%Y-%m-%d"))) %>% 
  group_by(CATCH_REGION, MGMT_AREA, FISHING.MONTH, jDay, FISHING.YEAR) %>% 
  summarize(catch = sum(CHINOOK_KEPT),
            boatDays = sum(VESSELS_OP)) %>% 
  select(catchReg = CATCH_REGION, area = MGMT_AREA, year = FISHING.YEAR, 
         month = FISHING.MONTH, jDay, catch, boatDays) %>%
  # filter(!boatDays == 0) %>% 
  ungroup() %>% 
  mutate(catchReg = as.character(catchReg),
         sumCPUE = catch / boatDays)

# write.csv(dailyCatch, here::here("data", "gsiCatchData", "commTroll",
#                                  "dailyCatch_WCVI.csv"), row.names = FALSE)

ggplot(data = dailyCatch %>% filter(!boatDays == 0)) +
  geom_point(aes(x = jDay, y = sumCPUE)) +
  facet_wrap(~catchReg) +
  samSim::theme_sleekX()


## Check anomalously large CPUE
# dailyCatch %>% 
#   filter(sumCPUE > 900)
# 
# areaCatch %>% 
#   filter(AREA_NAME == "management area 123", FISHING.YEAR == "2013", 
#          FISHING.MONTH == "5")
