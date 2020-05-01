## Map of PFMAs 
# April 30, 2020

library(tidyverse)
library(ggplot2)
library(maptools)
library(here)
library(rmapshaper)
library(geojsonio)
library(mapdata)


theme_set(ggsidekick::theme_sleek())

# relevant stat areas based on composition data 
areas <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                           "pstAggRollUpCatchProb.RDS")) %>% 
  pull(statArea) %>% 
  unique()

#big map
n_am <- map_data("worldHires", region = c("usa", "canada")) %>%
  filter(long < -110) %>%
  fortify(.)


#pfma polys 
pfma_simp_df <- readRDS(here("data", "gsiCatchData", "pfma", "trim_pfma_df.rds"))
# pfma_shp <- rgdal::readOGR(here("data", "gsiCatchData", "pfma"), 
#                             layer = "PFMA_Areas_50k") %>% 
#   sp::spTransform(., CRS("+proj=longlat +datum=WGS84"))
# pfma_simp <- ms_simplify(pfma_shp, sys = TRUE)
# pfma_simp_df <- pfma_simp %>% 
#   broom::tidy(pfma_simp, region = "STATAREA") %>% 
#   filter(lat < 51.5) %>% 
#   mutate(STATAREA = as.factor(id))
# saveRDS(pfma_simp_df, here("data", "gsiCatchData", "pfma", "trim_pfma_df.rds"))

# alternative version; ms_simplify is time consuming but worth it for reg.
#plot tweaks
# pfma <- broom::tidy(pfma_shp, region = "STATAREA")
# pfma_trim <- pfma %>% 
#   filter(lat < 51.5) %>% 
#   mutate(STATAREA = as.factor(id))

pfma_simp_df2 <- pfma_simp_df %>% 
  filter(STATAREA %in% areas) %>% 
  mutate(
    STATAREA = droplevels(STATAREA),
    stat_n = as.numeric(as.character(STATAREA)),
    catchReg = case_when(
      stat_n < 125 & stat_n > 27 ~ "SWVI",
      stat_n < 25 ~ "SWVI",
      TRUE ~ "NWVI"),
    statArea = fct_reorder(STATAREA, lat)
  )

#color palette
pal <- wesanderson::wes_palette("Rushmore1", 5, type = "continuous")
cols.named <- c(pal[1], pal[4])
cols.vals <- c('NWVI',  'SWVI')
alpha_labs <- levels(pfma_simp_df2$statArea)
alpha_vals <- rep(seq(0.7, 1, length = (length(alpha_labs) / 2)), 2)

ggplot() +
  coord_quickmap(xlim = c(-129, -122.2), ylim = c(47.75, 51), expand = TRUE,
                 clip = "on") +
  geom_polygon(data = pfma_simp_df2, aes(x = long, y = lat, group = group,
                                fill = catchReg, alpha = statArea),
               lwd = 1) +
  geom_polygon(data = n_am, aes(x=long, y = lat, group = group)) + 
  scale_fill_manual(name = "Region", labels = cols.vals, values = cols.named) +
  scale_alpha_manual(labels = alpha_labs, values = alpha_vals, guide = FALSE) +
  labs(x = "", y = "") +
  theme(plot.margin=unit(c(0.1,0,0,0), "mm"))

