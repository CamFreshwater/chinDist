## Make Chinook stock key based on high seas data
# July 5 2019

library(tidyverse)

adDat <- read.csv(here::here("data", "highSeas", "adultChinookHS.csv"))

stockKey <- adDat %>% 
  select(STOCK_1, REGION_1, HSS_REGION_1) %>% 
  distinct() 

write.csv(stockKey, here::here("data", "chinookStockKey.csv"), 
          row.names = FALSE)
