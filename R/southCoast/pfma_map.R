## Map of PFMAs 
# June 7, 2020

library(tidyverse)
library(ggplot2)
library(maptools)
library(here)
library(rmapshaper)
library(geojsonio)
library(mapdata)


theme_set(ggsidekick::theme_sleek())

# relevant stat areas based on composition data 
rec_areas <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(!is.na(area_n)) %>% 
  pull(area_n) %>% 
  unique()
# commercial data
comm_areas <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds")) %>% 
  pull(area_n) %>% 
  unique()
areas <- c(rec_areas, comm_areas)


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

# check to make sure all gsi areas are in pfma shape file after cropping
setdiff(areas, pfma_simp_df$STATAREA)

pfma_simp_df2 <- pfma_simp_df %>% 
  filter(STATAREA %in% areas) %>% 
  mutate(
    STATAREA = droplevels(STATAREA),
    stat_n = as.numeric(as.character(STATAREA)),
    region = case_when(
      stat_n > 124 ~ "NWVI",
      stat_n < 28 & stat_n > 24 ~ "NWVI",
      stat_n %in% c("20", "121", "21") ~ "Juan de Fuca",
      is.na(stat_n) ~ "Juan de Fuca",
      stat_n < 125 & stat_n > 120 ~ "SWVI",
      stat_n < 25 & stat_n > 20 ~ "SWVI",
      stat_n < 20 & stat_n > 12 ~ "Georgia Strait",
      stat_n %in% c("28", "29") ~ "Georgia Strait",
      stat_n %in% c("10", "11", "12", "111") ~ "Johnstone Strait",
      TRUE ~ NA_character_
    ),
    # catchReg = case_when(
    #   stat_n < 125 & stat_n > 27 ~ "SWVI",
    #   stat_n < 25 ~ "SWVI",
    #   TRUE ~ "NWVI"),
    region = fct_relevel(region, 'NWVI', 'SWVI', 'Johnstone Strait', 
                         'Juan de Fuca', 'Georgia Strait'),
    statArea = fct_reorder2(STATAREA, lat, region)
  )

pfma_simp_df2 %>% 
  select(region, statArea) %>% 
  distinct() %>% 
  arrange(statArea)


#color palette
# col.reg <- levels(pfma_simp_df2$region)
# pal <- wesanderson::wes_palette("BottleRocket2", 
#                                 length(col.reg), 
#                                 type = "continuous")
alpha_labs <- levels(pfma_simp_df2$statArea)
alpha_vals <- rep(seq(0.5, 0.8, length = (length(alpha_labs) / length(col.reg))),
                  length(col.reg))

ggplot() +
  coord_quickmap(xlim = c(-128.75, -122.2), ylim = c(48.25, 51), expand = TRUE,
                 clip = "on") +
  geom_polygon(data = pfma_simp_df2, aes(x = long, y = lat, group = group,
                                         colour = region, fill = region, 
                                         alpha = statArea),
               lwd = 1) +
  geom_polygon(data = n_am, aes(x=long, y = lat, group = group)) + 
  scale_colour_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  # scale_fill_manual(name = "Region", labels = col.reg, values = pal) +
  scale_alpha_manual(labels = alpha_labs, values = alpha_vals, guide = FALSE) +
  labs(x = "", y = "") +
  theme(plot.margin=unit(c(0.1,0,0,0), "mm"))

