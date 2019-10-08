## Stock-specific catch distributions from high seas survey
# Divide into juveniles (<250), immature (<450) and adult (>450)
# Restricted to stock ID'd samples
# August 2, 2019

library(tidyverse)
library(ggplot2)
library(ggmap)

source(here::here("R", "functions", "hotMapFunc.R"))

chinRaw <- read.csv(here::here("data", "highSeas", "fullChinookHS.csv"), 
                    stringsAsFactors = FALSE) %>% 
  filter(
    !PROB_1 < 0.5, #remove low probability GSI samples
  )

chin <- chinRaw %>%
  mutate(stationID = sapply(strsplit(FISH_NUMBER, "-"), function(x) {
    paste(x[1], x[2], x[3], sep = "-")
  })) %>%
  select(stationID, catchReg = REGION, lat = START_LAT, long = START_LONG,
         year = YEAR, month = MONTH, fishNumber = FISH_NUMBER, fl = SHIP_FL,
         stock = STOCK_1, reg = REGION_1, agg = HSS_REGION_1) %>% 
  mutate(month =  fct_relevel(as.factor(month), "FEB", "MAR", "APR", "MAY", 
                              "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"),
         agg = case_when(
           reg == "SOMN" ~ "BC SOUTH COAST",
           TRUE ~ agg),
         catchReg = fct_recode(catchReg, SoG = "GEORGIA STRAIT", 
                               QCSt = "QUEEN CHARLOTTE STRAIT", 
                               InBC = "INSIDE BC", JS = "JOHNSTONE STRAIT", 
                               QCSd = "QUEEN CHARLOTTE SOUND", 
                               WCVI = "VANCOUVER ISLAND", SEAK = "SE ALASKA", 
                               InSEAK = "INSIDE SE ALASKA", 
                               InVI = "INSIDE VANCOUVER ISLAND", 
                               DE = "DIXON ENTRANCE", 
                               SCAK = "SOUTH CENTRAL ALASKA",
                               HS = "HECATE STRAIT",  JDF = "JUAN DE FUCA", 
                               QCI = "QUEEN CHARLOTTE ISLANDS", 
                               PS = "PUGET SOUND", 
                               JDF_US = "JUAN DE FUCA - USA SIDE", 
                               DI = "DISCOVERY ISLANDS", 
                               BA = "BROUGHTON ARCHEPELIGO", 
                               KI = "KNIGHT INLET"))

# write.csv(chin, here::here("data", "highSeas", "cleanChinook.csv"),
#           row.names = FALSE)
  
## Visualize sampling effort across sampling regions
dropRegions <- chin %>% 
  group_by(catchReg) %>% 
  summarize(yearsFished = length(unique(year))) %>% 
  filter(yearsFished < 10)
dropReg <- dropRegions$catchReg

setCount <- chin %>% 
  filter(!catchReg %in% dropReg) %>% 
  mutate(season = case_when(
    month %in% c("MAY", "JUN", "JUL", "AUG") ~ "summer",
    month %in% c("SEP", "OCT", "NOV") ~ "fall",
    month %in% c("DEC", "JAN", "FEB", "MAR", "APR") ~ "winter")) %>% 
  group_by(year, season, catchReg) %>%
  tally(length(unique(stationID)), name = "setsWFish") %>% 
  ungroup()

pdf(here::here("figs", "highSeasCatch", "catchRegionSetCount.pdf"),
    height = 7, width = 8)
ggplot(setCount, aes(x = catchReg, y = year, fill = setsWFish)) + 
  geom_tile() + 
  viridis::scale_fill_viridis(name = "Sets", option = "viridis") +
  samSim::theme_sleekX() + 
  facet_wrap(~season, ncol = 1)
dev.off()


## Visualize annual catches for regional roll ups
summRegCatch <- chin %>% 
  filter(fl < 300) %>% 
  group_by(year, month, agg) %>% 
  tally() %>% 
  group_by(month, agg) %>% 
  summarize(mean = mean(n),
            median = median(n),
            sd = sd(n))

pdf(here::here("figs", "highSeasCatch", "catchCounts.pdf"),
    height = 7, width = 8)
ggplot(summRegCatch, aes(x = month, y = agg, fill = log(mean))) + 
  geom_tile() + 
  viridis::scale_fill_viridis(name = "Mean Catch", option = "viridis") +
  samSim::theme_sleekX()
dev.off()

## Generate heat maps
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
