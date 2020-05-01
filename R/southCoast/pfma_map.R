## Map of PFMAs 
# April 30, 2020

library(tidyverse)
library(ggplot2)
library(maptools)
library(here)
library(rmapshaper)
library(geojsonio)


theme_set(ggsidekick::theme_sleek())

#big map
n_am <- map_data("world", region = c("usa", "canada")) %>%
  filter(long < -110) %>%
  fortify(.)

#pfma polys 
pfma_shp <- rgdal::readOGR(here("data", "gsiCatchData", "pfma"), 
                            layer = "PFMA_Areas_50k") %>% 
  sp::spTransform(., CRS("+proj=longlat +datum=WGS84"))
pfma_simp <- ms_simplify(pfma_shp, sys = TRUE)
pfma_simp_df <- pfma_simp %>% 
  broom::tidy(pfma_simp, region = "STATAREA") %>% 
  filter(lat < 51.5) %>% 
  mutate(STATAREA = as.factor(id))
saveRDS(pfma_simp_df, here("data", "gsiCatchData", "pfma", "trim_pfma_df.rds"))

# alternative version; ms_simplify is time consuming but worth it for reg.
#plot tweaks
# pfma <- broom::tidy(pfma_shp, region = "STATAREA")
# pfma_trim <- pfma %>% 
#   filter(lat < 51.5) %>% 
#   mutate(STATAREA = as.factor(id))

ggplot() +
  geom_polygon(data = n_am, aes(x=long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-131, -122), ylim = c(47.5, 51.5), ratio = 1.3) +
  geom_polygon(data = pfma_simp_df, aes(x = long, y = lat, group = group,
                                fill = STATAREA), 
               lwd = 1, alpha = 0.4) 

world.plot <- ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
  coord_quickmap(xlim = c(30, 170), ylim = c(-25, 20), expand = TRUE,
                 clip = "on") + 
  labs(x = '', y = '') +
  geom_rect(data = bbox, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax=ymax, group=island), fill='transparent',lwd=0.5, col='red') +
  theme(plot.margin=unit(c(0.1,0,0,0), "mm"))


crit_hab <- rgdal::readOGR(here::here("data", "criticalHabitatShapeFiles"), 
                           layer = "Proposed_RKW_CriticalHabitat update_SWVI_CSAS2016") %>% 
  #transform to WGS84
  sp::spTransform(., sp::CRS("+proj=longlat +datum=WGS84"))