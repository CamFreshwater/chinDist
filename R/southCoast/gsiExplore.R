## Explore catch data from southcoast GSI
# August 22, 2019
# Goal is to identify stock-specific catches in NWVI and SWVI at monthly 
# and yearly time scales

library(tidyverse); library(ggplot2)

catch <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                             "stockCatch_WCVI.csv"))
