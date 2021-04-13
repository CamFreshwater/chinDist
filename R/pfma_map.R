## Map of PFMAs 
# June 7, 2020

library(ggplot2)
library(maptools)
library(here)
library(rmapshaper)
library(geojsonio)
library(mapdata)
library(raster)
library(maps) 
library(mapdata)
library(rgeos)
library(rgdal)
library(ggsn)
library(tidyverse)

theme_set(ggsidekick::theme_sleek())

# relevant stat areas based on composition data 
areas <- readRDS(here::here("data", "gsiCatchData", "pfma", "areas_to_plot.RDS"))

#pfma polys 
# pfma_area_simp_df <- readRDS(here("data", "gsiCatchData", "pfma", "trim_pfma_df.rds"))
# pfma_shp1 <- rgdal::readOGR(here("data", "gsiCatchData", "pfma", "shape_files"),
#                             layer = "PFMA_Areas_50k") %>%
#   sp::spTransform(., CRS("+proj=longlat +datum=WGS84"))
# pfma_simp <- ms_simplify(pfma_shp1, sys = TRUE)
# pfma_simp_df <- pfma_simp %>%
#   broom::tidy(pfma_simp, region = "STATAREA") %>%
#   filter(lat < 51.5) %>%
#   mutate(STATAREA = as.factor(id))
# saveRDS(pfma_simp_df, here("data", "gsiCatchData", "pfma", "trim_pfma_df.rds"))

#subarea polys
pfma_simp_df <- readRDS(here("data", "gsiCatchData", "pfma", 
                             "trim_subarea_df.rds"))
# pfma_shp <- rgdal::readOGR(here("data", "gsiCatchData", "pfma", "shape_files"),
#                            layer = "PFMA_Subareas_50k") %>% 
#   sp::spTransform(., CRS("+proj=longlat +datum=WGS84"))
# pfma_subarea_simp_df <- pfma_simp %>%
#   broom::tidy(pfma_simp, region = "LABEL") %>%
#   filter(lat < 51.5) %>%
#   mutate(subarea = as.factor(id))
# pfma_subarea_simp_df$pfma <- stringr::str_split(pfma_subarea_simp_df$subarea, "-") %>%
#   lapply(., function (x) x[1]) %>%
#   unlist()
# saveRDS(pfma_subarea_simp_df,
#         here("data", "gsiCatchData", "pfma", "trim_subarea_df.rds"))


# check to make sure all gsi areas are in pfma shape file after cropping
setdiff(areas, pfma_simp_df$pfma)

pfma_simp_df2 <- pfma_simp_df %>% 
  filter(pfma %in% areas) %>% 
  mutate(
    #adjust  subareas 13-1 and 13-2 to be area 14 (NSoG) consistent with changes
    #to model
    pfma = case_when(
      subarea %in% c("13-1", "13-2") ~ "14",
      TRUE ~ pfma
    ),
    stat_n = as.numeric(as.character(pfma)),
    region = case_when(
      stat_n > 124 ~ "NWVI",
      stat_n < 28 & stat_n > 24 ~ "NWVI",
      stat_n %in% c("20", "121", "21") ~ "Juan de Fuca Strait",
      is.na(stat_n) ~ "Juan de Fuca Strait",
      stat_n < 125 & stat_n > 120 ~ "SWVI",
      stat_n < 25 & stat_n > 20 ~ "SWVI",
      stat_n %in% c("14", "15", "16") ~ "N. Strait of Georgia",
      stat_n %in% c("17", "18", "19", "28", "29") ~ "S. Strait of Georgia",
      stat_n %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits",
      TRUE ~ NA_character_
    ),
    region = fct_relevel(region, 'NWVI', 'SWVI', 
                         'Queen Charlotte and\nJohnstone Straits',
                         'N. Strait of Georgia', 'S. Strait of Georgia',
                         'Juan de Fuca Strait'),
    statArea = fct_reorder2(as.factor(pfma), lat, region)
  )

pfma_simp_df2 %>% 
  select(region, statArea) %>% 
  distinct() %>% 
  arrange(region)

# high res shapefile of BC coastline
# coast <- map_data("worldHires", region = c("usa", "canada")) %>%
#   filter(long < -110) %>%
#   fortify(.)
# coast_shp <- readOGR(here("data", "coast_shapefiles", "COAST_TEST2.shp"))
# sBC_crop <- extent(-128.85, -122.1, 48.2, 51.2)
# coast2 <- crop(coast_shp, sBC_crop) %>% 
#   fortify(.)
# coast_shp_simp <- crop(coast_shp, sBC_crop) %>% 
#   ms_simplify(., sys = TRUE)
# coast2b <- fortify(coast_shp_simp)
# saveRDS(coast2b, here("data", "gsiCatchData", "pfma", "trim_coastline_df.rds"))
coast2b <- readRDS(here::here("data", "gsiCatchData", "pfma", 
                              "trim_coastline_df.rds"))


#color palette
col.reg <- levels(pfma_simp_df2$region)
brew_pal <- RColorBrewer::brewer.pal(n = length(col.reg), "Dark2")
names(brew_pal) <- col.reg
saveRDS(brew_pal, here::here("generated_data", "color_pal.RDS"))
disco_pal <- disco::disco(palette = "bright", n = length(col.reg))[c(2,1,3:6)]
names(disco_pal) <- col.reg
saveRDS(disco_pal, here::here("generated_data", "disco_color_pal.RDS"))

alpha_labs <- levels(pfma_simp_df2$statArea)
alpha_vals <- rep(seq(0.6, 1, length = (length(alpha_labs) / length(col.reg))),
                  length(col.reg))

pfma_map <- ggplot() +
  coord_quickmap(xlim = c(-128.75, -122.8), ylim = c(48.25, 51), expand = TRUE,
                 clip = "on") +
  geom_polygon(data = pfma_simp_df2, aes(x = long, y = lat, group = group,
                                         colour = region,
                                         fill = region,
                                         alpha = statArea),
               lwd = 0.1) +
  geom_polygon(data = coast2b, aes(x=long, y = lat, group = group)) + 
  scale_colour_manual(values = disco_pal, name = "", guide = FALSE) +
  scale_fill_manual(name = "Region", labels = col.reg, values = disco_pal) +
  scale_alpha_manual(labels = alpha_labs, values = alpha_vals, guide = FALSE) +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  theme(plot.margin=unit(c(0.1,0,0,0), "mm"),
        legend.position = "top")

png(here::here("figs", "ms_figs", "main", "pfma_map.png"), 
    res = 400, units = "in", height = 4, width = 6)
pfma_map
dev.off()

# pdf(here::here("figs", "ms_figs", "pfma_map.pdf"))
# pfma_map
# dev.off()

saveRDS(pfma_map, here::here("generated_data", "pfma_map.rds"))
