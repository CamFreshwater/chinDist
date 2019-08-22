## Explore catch data from southcoast GSI
# August 22, 2019
# Goal is to identify stock-specific catches in NWVI and SWVI at monthly 
# and yearly time scales

library(tidyverse); library(ggplot2)
source(here::here("R", "functions", "pooledVar.R"))

catch <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                             "stockCatch_WCVI.csv"))

regCatch <- catch %>% 
  group_by(Region2Name, catchReg, month, year, sumEffort, labN) %>% 
  summarize(regCatch = sum(estCatch),
            regCatchVar = sum(varCatch)) %>% 
  mutate(regCPUE = regCatch / sumEffort,
         varCPUE = regCatchVar * (sumEffort)^-2)

ggplot(regCatch, aes(x = as.factor(month), y = regCatch, fill = catchReg)) +
  geom_boxplot() +
  facet_wrap(~Region2Name, scales = "free_y")

ggplot(regCatch, aes(x = as.factor(month), y = regCPUE, fill = catchReg)) +
  geom_boxplot() +
  facet_wrap(~Region2Name, scales = "free_y")

xx <- rnorm(1000, 0, (3.56^2))

## Focus on BC stocks
bcStocks <- unique(catch$Region2Name)[c(1:3, 5:7)]
bcCatch <- catch %>%
  filter(Region2Name %in% bcStocks) %>% 
  group_by(Region1Name, catchReg, month, year, sumEffort, labN) %>% 
  summarize(regCatch = sum(estCatch),
            regCatchVar = sum(varCatch)) %>% 
  mutate(regCPUE = regCatch / sumEffort,
         varCPUE = regCatchVar * (sumEffort)^-2) %>% 
  ungroup() 

ggplot(bcCatch, aes(x = as.factor(month), y = regCPUE, fill = catchReg)) +
  geom_boxplot() +
  facet_wrap(~Region1Name, scales = "free_y")



## Bootstrap BC stocks to account for sampling uncertainty 
# TO BE COMPLETED
annualBCCatch <- bcCatch %>% 
  group_by(Region1Name, catchReg, month) %>% 
  summarize(meanRegCatch = mean(regCatch),
            pooledRegCatchVar = pooledVar(regCatchVar, labN))

bcList <- split(annualBCCatch[1:2,], seq(nrow(annualBCCatch[1:2,])))

simOut <- lapply(bcList, function(x) {
 data.frame(Region1Name = x$Region1Name,
            catchReg = x$catchReg,
            month = x$month,
            trial = seq(from = 1, to = 10000, by = 1),
            catch = rnorm(10000, x$meanRegCatch, x$pooledRegCatchVar^2))
})

