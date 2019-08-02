## Stock-specific catch distributions from high seas survey
# Divide into juveniles (<250), immature (<450) and adult (>450)
# Restricted to stock ID'd samples
# August 2, 2019

library(tidyverse)
library(ggplot2)
library(ggmap)

chinRaw <- read.csv(here::here("data", "highSeas", "fullChinookHS.csv"), 
                    stringsAsFactors = FALSE) %>% 
  filter(
    !PROB_1 < 0.5, #remove low probability GSI samples
  )

chin <- chinRaw %>%
  mutate(stationID = sapply(strsplit(FISH_NUMBER, "-"), function(x) {
    paste(x[1], x[2], x[3], sep = "-")
  })) %>%
  select(stationID, lat = START_LAT, long = START_LONG,
         year = YEAR, month = MONTH, fishNumber = FISH_NUMBER, fl = SHIP_FL,
         stock = STOCK_1, reg = REGION_1, agg = HSS_REGION_1) %>% 
  mutate(month =  fct_relevel(as.factor(month), "FEB", "MAR", "APR", "MAY", 
                              "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"))

hotMap <- function(fishDat, longRange = c(-129.5, -123), 
                   latRange = c(48, 52), facet = NULL) {
  nAm <- map_data("world") %>% 
    filter(region %in% c("Canada", "USA"))
  
  p <- ggplot(fishDat) +
    stat_density2d(aes(x = long, y = lat, fill = ..level..), 
                   geom = "polygon",
                   alpha = 0.5) +
    scale_fill_gradient(low = "green", high = "red") +
    scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
    lims(x = longRange, y = latRange) +
    geom_polygon(data = nAm, aes(x = long, y = lat, group = group), 
                 color = "black", fill = "gray80") 
  if (facet == "month") {
    p <- p +
      facet_wrap(~month)
  } 
  return(p)
}
  
hotMap(chin, longRange = c(-140, -120), latRange = c(47, 57))


upFr <- chin %>% 
  filter(reg == "UPFR")

hotMap(upFr, longRange = c(-140, -120), latRange = c(47, 57), facet = "month") 
