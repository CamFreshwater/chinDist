## Raw catch composition
# January 10, 2019
# Use raw assignment probabilities to generate catch composition estimates at
# regional scale. Precursor to fitting multinomial models and generating 
# estimates of predicted catch. 

library(tidyverse)
library(ggplot2)

reg3 <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                            "reg3RollUpCatchProb.csv"))
