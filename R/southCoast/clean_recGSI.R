## Explore rec GSI data
# Preliminary investigation into structure and patterns of recreational GSI
# data. Original data in chin directory of local (imported from salmon drive
# on May 15): southCoastDatabase/chinoo_rec_downloadApr23.xlsx

library(tidyverse)

rec_full <- read.csv(here::here("data", "gsiCatchData", "rec", 
                               "rec_gsi_may2020.csv"), stringsAsFactors = F)


# pull stocks to add to stockkey repo
# stk_out <- rec_full %>% 
#   filter(!DNA_RESULTS_STOCK_1 == "",
#          !is.na(PROB_1)) %>% 
#   select(s1 = DNA_RESULTS_STOCK_1, s2 = DNA_STOCK_2, s3 = DNA_STOCK_3, 
#          s4 = DNA_STOCK_4, s5 = DNA_STOCK_5, REGION_1_ROLLUP) %>% 
#   pivot_longer(., cols = s1:s5, names_to = "rank", values_to = "stock") %>%
#   select(stock, sc_reg1 = REGION_1_ROLLUP) %>% 
#   filter(!stock == "") %>% 
#   distinct()
# saveRDS(stk_out, here::here("data", "stockKeys", "rec_gsi_stocks.rds"))

stockKey <- readRDS(here::here("data", "stockKeys", "finalStockList_May2020.rds"))

# data frame of probabilities
temp_prob <- rec_dat %>% 
  filter(!DNA_RESULTS_STOCK_1 == "") %>% 
  rename(p1 = PROB_1, p2 = PROB_2, p3 = PROB_3, p4 = PROB_4, p5 = PROB_5) %>% 
  pivot_longer(., cols = c(p1, p2, p3, p4, p5),
               names_to = "rank_prob", values_to = "prob") %>%  
  select(BIOKEY, COLLECTION_DATE, rank_prob, prob)
  
# data frame of stock IDs
rec_long <- rec_dat %>% 
  filter(!DNA_RESULTS_STOCK_1 == "") %>% 
  rename(s1 = DNA_RESULTS_STOCK_1, s2 = DNA_STOCK_2, s3 = DNA_STOCK_3, 
         s4 = DNA_STOCK_4, s5 = DNA_STOCK_5,
         p1 = PROB_1, p2 = PROB_2, p3 = PROB_3, p4 = PROB_4, p5 = PROB_5) %>% 
  pivot_longer(., cols = c(s1, s2, s3, s4, s5), 
               names_to = "rank", values_to = "stock") %>% 
  cbind(., temp_prob %>% select(rank_prob, prob)) %>% 
  #add regional roll ups
  left_join(., stockKey, by = "stock") 


## subset and remove extra columns to match format of commercial data
# small sample size pfmas to drop
rec_gsi <- rec_long %>% 
  filter(!is.na(prob)) %>%
  mutate(date = as.Date(as.numeric(as.character(DAYOFYEAR - 1)),
                        origin = as.Date(paste(YEAR, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
         week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d"))
  ) %>% 
  #adjust probabilities 
  group_by(BIOKEY) %>% 
  mutate(total_prob = sum(prob),
         adj_prob = prob / total_prob) %>% 
  ungroup() %>% 
  select(area = PFMA, subarea = SUBAREA, year = YEAR, gear = SAMPLE_TYPE, 
         year_day = DAYOFYEAR, fish_num = FISH_NO, date, month, week, 
         id = BIOKEY, 
         fl = LENGTH_MM, sex = SEX, prob, stock, adj_prob, Region1Name:pst_agg)  


saveRDS(rec_gsi, here::here("data", "gsiCatchData", "commTroll",
                         "recIndProbsLong.rds"))


# random exploratory
rec_dat <- rec_full %>% 
  filter(!RESOLVED_STOCK_ORIGIN == "") %>% 
  mutate(MONTH = fct_reorder(MONTH, DAYOFYEAR))

# sample distribution by months
rec_dat %>% 
  ggplot(.) +
  geom_bar(aes(x = MONTH)) +
  facet_wrap(~PFMA)

# stock composition of samples
rec_dat %>% 
  group_by(RESOLVED_STOCK_ROLLUP, RESOLVED_STOCK_SOURCE) %>% 
  tally() %>% 
  print(n = Inf)

rec_dat %>% 
  filter(RESOLVED_STOCK_ROLLUP %in% c("Fraser Spring 4.2", 
                                      "Fraser Spring 5.2",
                                      "Fraser Summer 4.1")) %>% 
  group_by(RESOLVED_STOCK_ROLLUP, RESOLVED_STOCK_ORIGIN) %>% 
  tally() %>% 
  print(n = Inf)
