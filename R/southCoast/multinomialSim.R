### Simulated stock composition data
## Generate stock composition data with various structural properties and then
# recover with models 
## Oct. 4 2019

require(tidyverse)
require(rethinking)


## First simulate five year's composition with three stocks
yrs <- seq(1, 5, by = 1)
ppnStocks <- c(0.2, 0.4, 0.6)
nStocks <- length(ppnStocks)
error <- runif(nStocks, 0.0001, 0.9999)
samSim::ppnAgeErr(ppnStocks, tau, error)


N <- 500 #number of fish 

