## Explore southcoast GSI Access database
# Provided by Bryan Rusch Aug 19, 2019

library(RODBC)

# Access database saved locally to work computer
file_path <- "C:/Users/FRESHWATERC/Documents/chinook/southCoastDatabase/southCoastGSICatch.accdb"
con <- odbcConnectAccess2007(file_path)

# Check out table names
sqlTables(con, tableType = "TABLE")$TABLE_NAME

# Useful queries provided by Bryan
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

gsiQry <- "SELECT [Area G Commercial Stock Comp Control Table].[CHINOOK YEAR],
                  [Area G Commercial Stock Comp Control Table].[FISHING YEAR], 
                  [Area G Commercial Stock Comp Control Table].[FISHING MONTH], 
                  [Area G Commercial Stock Comp Control Table].CATCH_REGION, 
                  [Area G Catch by Year, Month and Region].SumOfCHINOOK_KEPT, 
                  [Region 3 Sample Results].[Region 3 Code], 
                  [Region 3 Sample Results].Estimate, 
                  [Region 3 Sample Results].SD, 
                  ([Estimate]/100)*[SumOFCHINOOK_KEPT] AS Catch, 
                  (([SD]/100)*[SumOfCHINOOK_KEPT])^2 AS CatchVar
           FROM [Area G Commercial Stock Comp Control Table]
           INNER JOIN [Region 3 Sample Results] 
           ON [Area G Commercial Stock Comp Control Table].[Sample ID] = 
              [Region 3 Sample Results].[Sample ID]) 
           WHERE ((([Area G Commercial Stock Comp Control Table].[Direct Sample])=Yes));
"
reg3 <- sqlQuery(con, gsiQry)
head(reg3)



# gsiQry <- "SELECT [Area G Commercial Stock Comp Control Table].[CHINOOK YEAR],
#                   [Area G Commercial Stock Comp Control Table].[FISHING YEAR], 
#                   [Area G Commercial Stock Comp Control Table].[FISHING MONTH], 
#                   [Area G Commercial Stock Comp Control Table].CATCH_REGION, 
#                   [Area G Catch by Year, Month and Region].SumOfCHINOOK_KEPT, 
#                   [Region 3 Sample Results].[Region 3 Code], 
#                   [Region 3 Sample Results].Estimate, 
#                   [Region 3 Sample Results].SD, 
#                   ([Estimate]/100)*[SumOFCHINOOK_KEPT] AS Catch, 
#                   (([SD]/100)*[SumOfCHINOOK_KEPT])^2 AS CatchVar
#            FROM [Area G Catch by Year, Month and Region] 
#             INNER JOIN ([Area G Commercial Stock Comp Control Table] 
#             INNER JOIN [Region 3 Sample Results] 
#             ON [Area G Commercial Stock Comp Control Table].[Sample ID] = 
#               [Region 3 Sample Results].[Sample ID]) 
#             ON ([Area G Catch by Year, Month and Region].CATCH_REGION = 
#               [Area G Commercial Stock Comp Control Table].CATCH_REGION) 
#             AND ([Area G Catch by Year, Month and Region].[FISHING MONTH] = 
#               [Area G Commercial Stock Comp Control Table].[FISHING MONTH]) 
#             AND ([Area G Catch by Year, Month and Region].[FISHING YEAR] = 
#               [Area G Commercial Stock Comp Control Table].[FISHING YEAR])
#            WHERE ((([Area G Commercial Stock Comp Control Table].[Direct Sample])=Yes));
# "