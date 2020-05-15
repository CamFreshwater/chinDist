## Explore rec GSI data
# Preliminary investigation into structure and patterns of recreational GSI
# data. Original data in chin directory of local (imported from salmon drive
# on May 15): southCoastDatabase/chinoo_rec_downloadApr23.xlsx

library(tidyverse)

rec_full <- read.csv(here::here("data", "gsiCatchData", "rec", 
                               "rec_gsi_may2020.csv"), stringsAsFactors = F)

rec_dat <- rec_full %>% 
  select(BIOKEY, ID2, YEAR:DAYOFYEAR, PFMA, SUBAREA, SIZE_CAT:RESOLVED_AGE,
         REGION_1_NAME:PROB_1) %>% 
  filter(!RESOLVED_STOCK_ORIGIN == "") %>% 
  mutate(MONTH = fct_reorder(MONTH, DAYOFYEAR))

# small sample size pfmas to drop
pfma_n <- rec_dat %>% 
  group_by(PFMA) %>% 
  tally()
pfma_keep <- pfma_n %>%
  filter(!n < 50) %>% 
  pull(PFMA)

rec_dat <- rec_dat %>% 
  filter(PFMA %in% pfma_keep)

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


# DATA CLEAN -------------------------------------------------------------------
# convert to format equivalent to commercial catch


