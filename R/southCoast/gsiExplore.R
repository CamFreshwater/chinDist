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
               FROM [Area G FOS Catch Estimates] INNER JOIN [Management Regions] 
                ON [Area G FOS Catch Estimates].MGMT_AREA = 
                  [Management Regions].[Managment Area]
               WHERE ((([Area G FOS Catch Estimates].TARGETS_CHINOOK)='Yes'));
"
areaCatch <- sqlQuery(con, catchQry)
head(areaCatch)

# all gsi samples
sampInvQry <- "SELECT *
               FROM [Sample Inventory]"
sampInv <- sqlQuery(con, sampInvQry) %>% 
  select(-F14, -F15, -F16) %>% 
  rename_all(list(~make.names(.))) %>% #get rid of spaces in col names
  filter(!Fishery == "Area G Subleg")
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
  left_join(., trimSampInv, by = "Sample.ID") 
# %>% 
  # select(Sample.ID:SD, Stock.Name:Lab.Reported.N)

dum <- stockComp %>% filter(is.na(Year))
dum71 <- stockComp %>% filter(Sample.ID == "71")

dum %>% 
  filter(is.na(Vial.ID)) %>% 
  select(Sample.ID) %>% 
  distinct()

regQry <- "SELECT *
           FROM [Region 2 Sample Results] 
           INNER JOIN [Region 2 Stock Names] 
           ON [Region 2 Sample Results].[Region 2 Code] = 
            [Region 2 Stock Names].Region2Code;"
regComp <- sqlQuery(con, regQry) %>% 
  rename_all(list(~make.names(.))) %>% 
  select(Sample.ID, Region2Code, Region2Name, Estimate, SD)
