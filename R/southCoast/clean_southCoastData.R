## Clean catch/effort data from southcoast 
# Commercial data provided by Bryan Rusch Aug 19, 2019 as Access database
# Note: the tables within the access database are not joined, limiting the
# utility of multitable queries; to ensure repeatability the majority of data
# cleaning occurs after tables have been imported to R.
# Recreational data provided by Wilf Luedke May 27, 2020 as Excel file
# Note: information on data selection for rec is in .txt file in directory

library(RODBC); library(tidyverse); library(ggplot2)


## CLEAN COMMERCIAL ------------------------------------------------------------

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
  filter(FISHING.YEAR %in% yrs) %>% 
  # correct anomalously low vessel operating after looking at FOS
  mutate(VESSELS_OP = 
           case_when(
             FISHING_DATE == "2013-05-15" & MGMT_AREA == "123" ~ 39,
             TRUE ~ VESSELS_OP)
         )
# write.csv(areaCatch, here::here("data", "gsiCatchData", "commTroll",
#                                 "fosCatch.csv"), 
#           row.names = FALSE)

#summary of how catch/effort data is distributed through time
# fosSumm <- areaCatch %>%
#   group_by(CATCH_REGION, FISHING.MONTH, FISHING.YEAR, MGMT_AREA) %>%
#   tally(name = "daysWithData")
# 
# write.csv(fosSumm, here::here("data", "gsiCatchData", "commTroll", 
#                               "fosSummary.csv"))

## Generate aggregate catch by Julian Day to match individual data from genetics
# lab
areaCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                                 "fosCatch.csv"))

dailyCatch <- areaCatch %>% 
  mutate(jDay = lubridate::yday(as.POSIXlt(FISHING_DATE, 
                                           format="%Y-%m-%d"))) %>% 
  group_by(CATCH_REGION, MGMT_AREA, FISHING.MONTH, jDay, FISHING.YEAR) %>% 
  #normally only KEPT chinook contribute to GSI samples, but certain
  #fisheries include live sampling where CHINOOK_RELD = catch
  summarize(catch =
              case_when(
                grepl("CN DNA Sampling", OPNG_DESC) ~ sum(CHINOOK_RELD),
                TRUE ~ sum(CHINOOK_KEPT)),
            boatDays = sum(VESSELS_OP)) %>% 
  rename(catchReg = CATCH_REGION, area = MGMT_AREA, year = FISHING.YEAR, 
         month = FISHING.MONTH) %>%
  ungroup() %>% 
  mutate(catchReg = as.character(catchReg),
         cpue = catch / boatDays) %>% 
  arrange(area, year, month, jDay)


saveRDS(dailyCatch, here::here("data", "gsiCatchData", "commTroll",
                                 "dailyCatch_WCVI.rds"))

ggplot(data = dailyCatch %>% filter(!boatDays == 0)) +
  geom_point(aes(x = jDay, y = sumCPUE)) +
  facet_wrap(~catchReg) +
  ggsidekick::theme_sleek()


## Check anomalously large CPUE
# dailyCatch %>% 
#   filter(sumCPUE > 900)
# 
# areaCatch %>% 
#   filter(AREA_NAME == "management area 123", FISHING.YEAR == "2013", 
#          FISHING.MONTH == "5")


# CLEAN REC DATA ---------------------------------------------------------------

rec_catch <- read.csv(here::here("data", "gsiCatchData", "rec",
                                 "southcoast_rec_cpue_Feb2020.csv")) %>% 
  mutate(strata = paste(PFMA, YEAR, MONTH, DISPOSITION, sep = "_"))

#initial explore
# spatio-temporal data coverage
rec_catch %>% 
  select(PFMA, YEAR, MONTH, DISPOSITION) %>%
  distinct() %>% 
  ggplot() + 
  geom_bar(aes(x = MONTH, fill = DISPOSITION)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~PFMA)

# estimate type breakdown
rec_catch %>% 
  group_by(STATUS) %>% 
  tally()

unpub_strata <- rec_catch %>% 
  filter(STATUS != "Published Estimate - Full Month") %>% 
  pull(strata)

rec_catch %>% 
  filter(STATUS == "Published Estimate - Full Month",
         strata %in% unpub_strata)
#none seem to be cross referenced so retain both

# clean
rec_catch1 <- rec_catch %>% 
  mutate(
    pfma_n = as.numeric(str_remove_all(PFMA, "PFMA ")),
    region = case_when(
      pfma_n > 124 ~ "NWVI",
      pfma_n < 28 & pfma_n > 24 ~ "NWVI",
      pfma_n %in% c("20", "121", "21") ~ "Juan de Fuca",
      is.na(pfma_n) ~ "Juan de Fuca",
      pfma_n < 125 & pfma_n > 120 ~ "SWVI",
      pfma_n < 25 & pfma_n > 20 ~ "SWVI",
      pfma_n < 20 & pfma_n > 12 ~ "Georgia Strait",
      pfma_n %in% c("28", "29") ~ "Georgia Strait",
      pfma_n %in% c("10", "11", "12", "111") ~ "Johnstone Strait"
    ),
    strata = paste(pfma_n, CREEL_SUB_AREA, YEAR, MONTH, sep = "_"),
    month = fct_relevel(MONTH, "January", "February", "March", "April", 
                        "May", "June", "July", "August", "September", 
                        "October", "November", "December"),
    month_n = as.numeric(month)
  ) %>%
  select(strata, month, month_n, year = YEAR, area = pfma_n,
         subarea = CREEL_SUB_AREA, region, 
         DISPOSITION, ADIPOSE_MARK, ESTIMATE, STANDARD_ERROR)

# separate effort and catch data, reformat, then recombine  
rec_eff <- rec_catch1 %>% 
  filter(DISPOSITION == "Effort") %>% 
  select(strata:region, mu_boat_trips = ESTIMATE, 
         se_boat_trips = STANDARD_ERROR)

rec_catch2 <- rec_catch1 %>% 
  filter(!DISPOSITION == "Effort") %>% 
  mutate(
    kept = case_when(
      DISPOSITION == "Kept" ~ "y",
      grepl("Released", DISPOSITION) ~ "n"
    ),
    legal = case_when(
      DISPOSITION == "Released Sub-Legal" ~ "sublegal",
      DISPOSITION %in% c("Released Legal", "Kept") ~ "legal"
      ),
    kept_legal = paste(kept, legal, sep = "_"),
    adipose_clip = case_when(
      ADIPOSE_MARK == "Adipose Marked" ~ "y",
      ADIPOSE_MARK == "Not Adipose Marked" ~ "n",
      TRUE ~ NA_character_)
    ) %>% 
  select(strata:region, legal, kept_legal, adipose_clip, mu_catch = ESTIMATE,
         se_catch = STANDARD_ERROR) %>% 
  distinct() %>% 
  filter(!is.na(mu_catch))

# rec_catch2 %>% 
#   group_by(strata, kept_legal, adipose_clip) %>% 
#   filter(n()>1)

# save list that includes uncertainty in both variables
rec_list <- list("catch" = rec_catch2, "effort" = rec_eff)

# combine into single dataframe that excludes uncertainty
rec_dat_out <- rec_catch2 %>% 
  select(-se_catch) %>%
  left_join(.,
            rec_eff %>% 
              select(strata, mu_boat_trips),
            by = "strata") %>% 
  select(-strata)

saveRDS(rec_list, here::here("data", "gsiCatchData", "rec",
                             "monthlyCatchList_rec.RDS"))
saveRDS(rec_dat_out, here::here("data", "gsiCatchData", "rec",
                             "monthlyCatch_rec.RDS"))
