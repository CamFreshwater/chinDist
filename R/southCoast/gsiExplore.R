## Explore southcoast GSI Access database
# Provided by Bryan Rusch Aug 19, 2019

library(RODBC); library(tidyverse)

# Access database saved locally to work computer
file_path <- "C:/Users/FRESHWATERC/Documents/chinook/southCoastDatabase/southCoastGSICatch.accdb"
con <- odbcConnectAccess2007(file_path)

# Check out table names
sqlTables(con, tableType = "TABLE")$TABLE_NAME

# Import majority of data directly from tables because they aren't linked with
# one another, limiting the utility of SQL queries that include joins

# all gsi samples
sampInvQry <- "SELECT *
               FROM [Sample Inventory]"
sampInv <- sqlQuery(con, sampInvQry) %>% 
  select(-F14, -F15, -F16) %>% 
  rename_all(list(~make.names(.))) #get rid of spaces in col names
trimSampInv <- sampInv %>% 
  select(Sample.ID, Vial.ID, Year, Fishery, Catch.Region, Catch.Area, 
         Month.Name, Lab.Reported.N)

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
  filter(Fishery == "Area G Com")
yrs <- unique(stockComp$Year)

# confirm one sample collected per catch region, year, and month
# dum <- stockComp %>% 
#   select(-c(Stock.Code:Vial.ID)) %>% 
#   distinct() %>% 
#   group_by(Catch.Region, Year, Month.Name) %>% 
#   tally()

regQry <- "SELECT *
           FROM [Stocks with Region Codes]
           LEFT JOIN [Region 2 Stock Names]
           ON [Stocks with Region Codes].Region2Code = 
              [Region 2 Stock Names].Region2Code;"
regKey <- sqlQuery(con, regQry) %>% 
  select(ID:FraserGroupCode, Region2Name)
head(regKey)

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
  filter(FISHING.YEAR %in% yrs)

#summarize catches by area, month and year to match GSI data
sumCatch <- areaCatch %>% 
  group_by(MGMT_AREA, FISHING.MONTH, FISHING.YEAR) %>% 
  summarize(sumCatch = sum(CHINOOK_KEPT),
            sumEffort = sum(VESSELS_OP)) %>% 
  select(MGMT_AREA, FISHING.MONTH, FISHING.YEAR, sumCatch, sumEffort)
head(sumCatch)

