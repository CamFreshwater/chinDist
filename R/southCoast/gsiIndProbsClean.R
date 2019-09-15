## Clean individual genetic data
# Data from west coast Vancouver Island commercial troll fishery
# Sep 14, 2019

library(tidyverse)
library(ggplot2)

dat <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbs_trim.csv"), 
                stringsAsFactors = FALSE)

# pull ID column, split and rename
xy <- data.frame(t(matrix(unlist(strsplit(as.vector(dat$sampleID), 
                                          split = " ")),
                        nrow = 5, ncol = length(dat$sampleID)))) %>% 
  rename("statArea" = X1, "year" = X2, "gear" = X3, "jDay" = X4, 
         "fishNum" = X5) 


dum <- sapply(strsplit(as.character(xy$year), "()"), function(x) print(x)[2:3])
