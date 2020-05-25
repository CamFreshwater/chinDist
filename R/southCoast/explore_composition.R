## Raw catch composition
# May 25, 2020
# Use raw assignment probabilities to generate catch composition estimates for both
# rec and commercial fisheries at various scales

library(tidyverse)

rec_gsi <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                              "recIndProbsLong.rds"))
comm_gsi <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                               "wcviIndProbsLong.rds"))

#combine
rec_trim <- rec_gsi %>% 
  mutate(fishery = "rec") %>% 
  select(id, area, year, month, week, adj_prob, stock, Region1Name, pst_agg, fishery)
gsi <- comm_gsi %>% 
  mutate(area = as.numeric(statArea),
         year = as.numeric(year),
         fishery = "comm") %>% 
  select(id, area, year, month, week, adj_prob = adjProb, stock, Region1Name, 
         pst_agg, fishery) %>% 
  rbind(., rec_trim) %>% 
  mutate(
    region = case_when(
      area > 124 ~ "NWVI",
      area < 28 & area > 24 ~ "NWVI",
      area < 125 & area > 120 ~ "SWVI",
      area < 25 & area > 20 ~ "SWVI",
      area < 20 & area > 12 ~ "Georgia Strait",
      area %in% c("28", "29") ~ "Georgia Strait",
      area %in% c("10", "11", "12", "111") ~ "Johnstone Strait",
      area %in% c("20") ~ "Juan de Fuca",
      is.na(area) ~ "Juan de Fuca"),
    reg1_trim = case_when(
      grepl("FR", pst_agg) ~ Region1Name,
      TRUE ~ "Non-Fraser"
    )
  )

#restrict to maximum assignments with cutoff based on PST aggregates
# helper function
calc_max_prob <- function(dat_in, thresh = 0.75) {
  out <- dat_in %>% 
    summarize(agg_prob = sum(adj_prob)) %>% 
    arrange(id, desc(agg_prob)) %>%
    ungroup() %>% 
    group_by(id) %>% 
    mutate(max_assignment = max(agg_prob)) %>% 
    # Remove samples where top stock ID is less than threshold probability
    filter(!agg_prob < max_assignment, 
           !max_assignment < thresh) %>% 
    ungroup() %>% 
    distinct() %>% 
    left_join(.,
              gsi %>% 
                select(id:week, fishery, region) %>% 
                distinct(),
              by = "id") %>% 
    select(-agg_prob, -max_assignment)
  
  colnames(out)[2] <- "agg"
  
  return(out)
}


# INITIAL EXPLORE --------------------------------------------------------------

pst_gsi <- gsi %>% 
  group_by(id, pst_agg) %>% 
  calc_max_prob(.)

reg1_gsi <- gsi %>% 
  group_by(id, reg1_trim) %>% 
  calc_max_prob(.)

# number of samples by area
pst_gsi %>% 
  # filter(!is.na(area)) %>% 
  ggplot(.) +
  geom_bar(aes(x = as.factor(month), fill = fishery)) +
  facet_wrap(~region)

# helper function to plot composition data based on aggregation
plot_comp <- function(gsi_in) {
  dum <- gsi_in %>%
    mutate(abb_reg = paste(abbreviate(region, 4), fishery, sep = "_"),
           strata = paste(abb_reg, month, sep = "_")) %>% 
    group_by(strata) %>% 
    mutate(nSampled = length(unique(id))) %>% 
    group_by(strata, agg) %>% 
    mutate(nByReg = length(unique(id)), 
           catch_prop = nByReg / nSampled) %>% 
    ungroup() %>%  
    select(strata, region, month, agg, catch_prop, nByReg, 
           nSampled, fishery, abb_reg) %>%
    distinct() 
  
  ggplot(dum) +
    geom_bar(aes(x = as.factor(month), y = catch_prop, fill = agg),
             stat = "identity") +
    scale_fill_viridis_d(option = "D") +
    labs(x = "", y = "Composition", fill = "Region of Origin") +
    facet_wrap(~abb_reg) + 
    ggsidekick::theme_sleek()
}

plot_comp(pst_gsi)
plot_comp(reg1_gsi)
