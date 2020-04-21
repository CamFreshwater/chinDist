## Clean CWT data
# Import and clean CWT data from RMIS and generate equivalent data to GSI ind
# probs clean
# Raw RMIS data query (should be equivalent to GSI) as follows:
# WHERE tag_status = '1'
# AND   species = '1'
# AND   recovery_date_year IN (2007,2008,2009,2010,2011,2012,2013,2014,2015)
# AND   fishery IN ('10','11','12','13','14','15','16','60')
# AND   recovery_location_rmis_region = 'WCVI'
# April 15, 2020

library(tidyverse)

## ADD STOCKS TO RECOVERY DATA -------------------------------------------------

datRaw <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                              "cwt_recov.txt"), 
                   stringsAsFactors = FALSE) %>% 
  filter(!grepl("50105", tag_code),
         !grepl("50104", tag_code),
         !grepl("60102", tag_code)) %>% 
  mutate(recovery_loc_code = str_squish(recovery_location_code))

# export tag codes to query rmis
tag_codes <- unique(stringr::str_pad(datRaw$tag_code, 6, pad = "0"))
write.table(tag_codes, here::here("data", "gsiCatchData", "commTroll", 
                                        "cwt_tag_code.txt"), 
            row.names = FALSE, col.names = FALSE)

# import generated tag_code key (from rmis)
key <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                           "cwt_tag_key.txt"), 
                stringsAsFactors = FALSE) %>% 
  mutate(
    stock_location_name = case_when(
      is.na(stock_location_name) ~ hatchery_location_name,
      TRUE ~ stock_location_name
    )
  )

#pull unique stock, region, basin combinations to combine with gsi in stockkey
#repo
key_out <- key %>%
  select(release = release_location_name,
         stock = stock_location_name,
         rmis_region = release_location_rmis_region,
         basin = release_location_rmis_basin,
         state = release_location_state) %>%
  distinct()
saveRDS(key_out, here::here("generatedData", "cwt_stock_key.RDS"))

# import completed stock_list
key_in <- readRDS(here::here("data", "stockKeys", 
                             "finalStockList_Apr2020.rds")) %>% 
  filter(id_type == "cwt") 

# add regional aggregates to key
key_final <- key %>% 
  select(tag_code = tag_code_or_release_id, 
         stock = stock_location_name) %>% 
  mutate(stock = toupper(stock)) %>% 
  left_join(., key_in, 
            by = c("stock"
                   )) 
  
# import recovery location key
recov_key <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                                 "cwt_recovery_location_key.csv"), 
                      stringsAsFactors = F) %>% 
  select(-Recoveries, -X, -X.1) %>% 
  mutate(recovery_loc_code = Code)

#merge stock ID and recovery location keys to cwt recoveries
rec_raw <- datRaw %>% 
  left_join(., key_final, by = "tag_code") %>% 
  left_join(., recov_key, by = "recovery_loc_code")


## FORMAT CWT DATA TO MATCH GSI ------------------------------------------------

rec_long <- rec_raw %>%
  mutate(statArea = NA,
         gear = "troll",
         date = as.POSIXct(as.character(rec_raw$recovery_date), 
                           format="%Y%m%d"),
         month = lubridate::month(date),
         week = lubridate::week(date), 
         jDay = lubridate::yday(date),
         season_c = case_when(
           month %in% c("12", "1", "2") ~ "w",
           month %in% c("3", "4", "5") ~ "sp",
           month %in% c("6", "7", "8") ~ "su",
           month %in% c("9", "10", "11") ~ "f"
         ),
         season = fct_relevel(season_c, "sp", "su", "f", "w"),
         month_n = as.numeric(month),
         month = as.factor(month_n),
         year =  as.factor(run_year),
         pres = 1
         ) %>% 
  select(year, month, week, jDay, gear, fishNum = recovery_id, date, 
         season_c, season, month_n,
         pres, catchReg = Basin, pst_agg)


#import preliminary gsi dataset for comparison

reg3 <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                         "reg3RollUpCatchProb.RDS"))
