## Clean CWT data
# Import and clean CWT data from RMIS and generate equivalent data to GSI ind
# probs clean
# April 15, 2020

library(tidyverse)

datRaw <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                              "cwt_recov.txt"), 
                   stringsAsFactors = FALSE) %>% 
  filter(!grepl("50105", tag_code),
         !grepl("50104", tag_code),
         !grepl("60102", tag_code))

# export tag codes to query rmis
tag_codes <- unique(stringr::str_pad(datRaw$tag_code, 6, pad = "0"))
write.table(tag_codes, here::here("data", "gsiCatchData", "commTroll", 
                                        "cwt_tag_code.txt"), 
            row.names = FALSE, col.names = FALSE)
unique(datRaw$tag_code)

# import generated tag_code key
key <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                           "cwt_tag_key.txt"), 
                stringsAsFactors = FALSE)

#pull unique stock, region, basin combinations to combine with gsi in stock key
key_out <- key %>% 
  select(release = release_location_name, hatchery = hatchery_location_name,
         stock = stock_location_name, rmis_region = release_location_rmis_region,
         basin = release_location_rmis_basin, state = release_location_state) %>% 
  distinct()
saveRDS(key_out, here::here("generatedData", "cwt_stock_key.RDS"))
