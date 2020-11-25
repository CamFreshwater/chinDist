## Combined model fits
# Nov 23, 2020
# Explore potential impacts of area 19 and 20 management measures.
# 1) Show contribution of different management action categories to samples by
# year, month and region.
# 2) Break down catch by kept/release.

library(tidyverse)

## COMPOSTION ------------------------------------------------------------------

# recreational composition data
rec <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(legal == "legal",
         region %in% c("Juan de Fuca Strait", "S. Strait of Georgia")) %>% 
  select(id, region, area, subarea, year, month = month_n, year_day = jDay, week,
         fl, 
         ad_clip, legal, release, sampler) %>% 
  distinct() %>% 
  droplevels() %>% 
  mutate(year_n = as.numeric(as.character(year)),
         fl = fl / 10
  ) %>% 
  #add restrictions based on size + clipping + spatio-temporal closures 
  #described in figure b4 in Dobson et al 2020
  mutate(
    restriction_group = case_when(
      subarea %in% c("19C", "19D", "19E") ~ "a",
      subarea %in% c("20C", "20D") ~ "b",
      subarea == "29DE" ~ "c",
      subarea == "20B" ~ "d",
      TRUE ~ NA_character_
    ),
    restriction_size = ifelse(fl > 67, "y", "n"),
    restriction_time = case_when(
      restriction_group == "a" & month > 2 & week < 24 ~ "y",
      restriction_group == "b" & month > 2 & week < 24  & year_n < 2018 ~ "y",
      restriction_group == "c" & month > 4 & week < 28 &  year_n < 2016 ~ "y",
      restriction_group == "c" & month > 2 & week < 28 & year_n < 2019 ~ "y",
      restriction_group == "c" & month > 2 & week < 32 & year_n < 2020 ~ "y",
      restriction_group %in% c("b", "d") & month > 6 & week < 28 & 
        year_n >= 2018 ~ "y",
      TRUE ~ NA_character_
    ),
    restricted = case_when(
      restriction_group == "a" & restriction_time == "y" &  
        restriction_size == "y" & ad_clip == "N" ~ "y",
      restriction_group == "b" & restriction_time == "y" &  
        restriction_size == "y" & ad_clip == "N" ~ "y",
      restriction_group == "c" & restriction_time == "y" ~ "y",
      # to account for more extreme restrictions in last two years
      restriction_group %in% c("b", "d") & restriction_time == "y" & 
        year_n >= 2018 ~ "y",
      TRUE ~ "n"
    ),
    man_category = case_when(
      is.na(restriction_group) ~ "status quo",
      !is.na(restriction_group) & restricted == "n" ~ "exempted",
      restricted == "y" ~ "restricted"
    )
  )


# stacked bar plot of samples available by management category
gsi_by_category <- rec %>% 
  filter(month > 2, month < 9) %>% 
  group_by(year, month, region, man_category) %>% 
  tally(name = "count") %>% 
  ggplot(., aes(x = as.factor(month), y = count, 
                fill = fct_reorder(man_category, -count))) + 
  geom_bar(position = "fill", stat  = "identity") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_d() +
  labs(x = "Month", y = "Proportion of Genetic Samples", 
       fill = "Management\nCategory") +
  facet_grid(year ~ region)

png(here::here("figs", "ms_figs", "comp_man_category.png"), res = 400, 
    units = "in",
    height = 5.5, width = 7)
gsi_by_category
dev.off()


## CATCH -----------------------------------------------------------------------

#recreational catch data - use subarea-month-year catch estimate to explicitly
#include release category
rec_catch <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "month_subarea_recCatch.RDS")) %>% 
  #drop sublegal fish and regions without genetics
  filter(legal == "legal",
         region %in% c("Juan de Fuca Strait", "S. Strait of Georgia")) %>% 
  #group by subarea to get rid of adipose and released legal duplicates
  group_by(month, month_n, year, area, subarea, region, kept_legal) %>% 
  summarize(subarea_catch = sum(mu_catch, na.rm = T),
            subarea_eff = mean(mu_boat_trips, na.rm = T),
            .groups = "drop") %>% 
  filter(!is.na(subarea_eff)) %>% 
  mutate(restriction_group = case_when(
    subarea %in% c("19C", "19D", "19E") ~ "a",
    subarea %in% c("20C", "20D") ~ "b",
    subarea == "29DE" ~ "c",
    subarea == "20B" ~ "d",
    TRUE ~ NA_character_
  ),
  restriction_time = case_when(
    restriction_group == "a" & month_n > 2 & month_n < 7 ~ "y",
    restriction_group == "b" & month_n > 2 & month_n < 7 & year < 2018 ~ "y",
    restriction_group == "c" & month_n > 4 & month_n < 8 & year < 2016 ~ "y",
    restriction_group == "c" & month_n > 2 & month_n < 8 & year < 2019 ~ "y",
    restriction_group == "c" & month_n > 2 & month_n < 9 & year < 2020 ~ "y",
    restriction_group %in% c("b", "d") & month_n > 6 & month_n < 7 & 
      year >= 2018 ~ "y",
    TRUE ~ NA_character_
  ),
  man_category = case_when(
    is.na(restriction_group) ~ "status quo",
    kept_legal == "y_legal" ~ "retained",
    kept_legal == "n_legal" ~ "released"
  )
  )

# bar plot of samples available
catch_by_category <- rec_catch %>% 
  filter(month_n > 2, month_n < 9,
         year > 2013) %>% 
  group_by(year, month, region, man_category) %>% 
  summarize(catch = sum(subarea_catch)) %>% 
  ggplot(., aes(x = month, y = catch, 
                fill = fct_reorder(man_category, -catch))) + 
  geom_bar(position = "fill", stat  = "identity") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_d(option = "D") +
  labs(x = "Month", y = "Proportion of Catch", 
       fill = "Management\nCategory") +
  facet_grid(year ~ region)

png(here::here("figs", "ms_figs", "catch_man_category.png"), res = 400, 
    units = "in",
    height = 5.5, width = 7)
catch_by_category
dev.off()
