# Quick breakdown of sampled BC stocks
# Generated for discussion with SARA biologists 
# Oct 28, 2020

library(tidyverse)

rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(legal == "legal") %>% 
  select(id, region, month, year, gear, adj_prob, stock, Region1Name, 
         Region3Name)

# commercial data
comm <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds")) %>% 
  # drop month-region-years where catch data are missing or where samples came 
  # from Taaq fishery, and no comm fishery was active in the same region
  mutate(temp_strata2 = paste(month, region, year, sep = "_")) %>% 
  select(id, region, month, year, gear, adj_prob, stock, Region1Name, 
         Region3Name)

key_stocks <- rbind(rec, comm) %>% 
  filter(grepl("OKA", stock) | Region3Name == "Fraser River") 

key_stocks %>% 
  group_by(Region1Name) %>% 
  summarize(inds = sum(adj_prob))
  
key_stocks %>% 
  filter(Region1Name == "Fraser_Fall") %>% 
  pull(stock) %>% 
  unique()
