## Stock-specific catch distributions from high seas survey
# Explore stock specific catches for juvenile Chinook during winter on WCVI
# Used to inform study design for potential sampling with south coast
# Nov. 4, 2019

library(tidyverse)
library(ggplot2)

chin <- read.csv(here::here("data", "highSeas", "cleanChinook.csv"))

## Focus on IVI Nov-March
winterVI <- chin %>% 
  filter(month %in% c("NOV", "DEC", "FEB", "MAR"),
         catchReg %in% c("WCVI", "InVI"))

ggplot(winterVI) +
  geom_boxplot(aes(x = as.factor(month), y = fl))

winterVI %>% 
  group_by(month, year) %>% 
  tally()

# Quick map
nAm <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))

ggplot(winterVI) +
  geom_point(aes(x = long, y = lat, color = reg)) +
  geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
           color = "black", fill = "gray80") + 
  coord_fixed(xlim = c(-129.5, -123), ylim = c(48, 52), ratio = 1.3) + 
  facet_wrap(~month)

# WCVI stocks and inlets only 
wcvi <- winterVI %>% 
  filter(agg == "WCVI",
         catchReg == "InVI")

ggplot(wcvi) +
  geom_bar(aes(x = stock)) +
  facet_wrap(~month)

# Quick map
ggplot(wcvi) +
  geom_point(aes(x = long, y = lat, color = stock)) +
  geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
           color = "black", fill = "gray80") + 
  coord_fixed(xlim = c(-128.25, -125), ylim = c(48.5, 51.6), ratio = 1.3) + 
  facet_wrap(~month)

wcvi %>% 
  group_by(stock) %>% 
  tally() %>% 
  print(n = nrow(.))
