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
    geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
             color = "black", fill = "gray80") + 
    lims(x = longRange, y = latRange) +
    samSim::theme_sleekX()
  if (!is.null(facet)) {
    if (facet == "month") {
      nSize <- fishDat %>% 
        group_by(month) %>% 
        tally()
      p <- p +
        facet_wrap(~month) 
    }
    if (facet == "year") {
      nSize <- fishDat %>% 
        group_by(year) %>% 
        tally()
      p <- p +
        facet_wrap(~year) 
    }
    p <- p +
      geom_text(data = nSize, 
                aes(x = -Inf, y = min(latRange), label = n), 
                vjust = 0, hjust = -1)
  }
  return(p)
}
  
hotMap(chin, longRange = c(-138, -120), latRange = c(47, 58.5))

## Iterate across regions of interest 
regions <- c("UPFR", "MUFR", "LWFR-F", "SOTH", "NOTH", "LWTH", "SOMN", "NOMN",
             "ECVI")
for (x in seq_along(regions)) {
  datJuv <- chin %>% 
    filter(reg == regions[x],
           fl < 250)
  pp <- hotMap(datJuv, longRange = c(-138, -120), latRange = c(47, 58.5), 
         facet = "month") +
    ggtitle(regions[x])
  pdf(here::here("figs", "highSeasCatch", "stockSpecific", "juv", 
                 paste(regions[x], "Juv.pdf", sep = "")),
      height = 7, width = 8)
  print(pp)
  dev.off()
  
  datAd <- chin %>% 
    filter(reg == regions[x],
           fl > 250)
  pp2 <- hotMap(datAd, longRange = c(-138, -120), latRange = c(47, 58.5), 
               facet = "month") +
    ggtitle(regions[x])
  pdf(here::here("figs", "highSeasCatch", "stockSpecific", "adult", 
                 paste(regions[x], "Adult.pdf", sep = "")),
      height = 7, width = 8)
  print(pp2)
  dev.off()
}


## Focus on CR spring
chin %>% 
  filter(agg == "COLUMBIA-SNAKE") %>% 
  select(reg) %>% 
  distinct()

juvCR <- chin %>% 
  # focus on spring/summer run
  filter(reg %in% c("MID COL-SP", "UP COL-SP", "UPPER COLUMBIA-SP", 
                    "UPPER COLUMBIA-SU/F", "MID COLUMBIA-SP", "SNAKE-SP/SU",
                    "UPPER WILLAMETTE", "UP WILLAMETTE"),
         fl > 150, fl < 400)

hotMap(juvCR, longRange = c(-132, -124.5), latRange = c(48, 52))
