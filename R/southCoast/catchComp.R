## Raw catch composition
# January 10, 2019
# Use raw assignment probabilities to generate catch composition estimates at
# regional scale. Precursor to fitting multinomial models and generating 
# estimates of predicted catch. 

library(tidyverse)
library(ggplot2)

reg3 <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                            "reg3RollUpCatchProb.RDS"))
fullGSI <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                         "wcviIndProbsLong_CLEAN.rds"))
weekCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                                  "weeklyCatch_WCVI.csv"), 
                      stringsAsFactors = TRUE)


## Calculate regional proportion assuming perfect assignment for each stat area-
#week and lumping all years; focus only 123/124 for now
max_prob <- reg3 %>% 
  group_by(flatFileID) %>% 
  mutate(max_assignment = max(aggProb),
         strata = paste(statArea, week, sep = "_")) %>% 
  # Remove samples where top stock ID is less than 75% probability
  filter(!aggProb < max_assignment, 
         !max_assignment < 0.75,
         statArea %in% c("123", "124"), 
         month %in% c(5, 6, 7, 8)) %>% 
  ungroup() %>% 
  distinct()

weekly_comp <- max_prob %>%
  group_by(strata) %>% 
  mutate(nSampled = length(unique(flatFileID))) %>% 
  group_by(strata, regName) %>% 
  mutate(nByReg = length(unique(flatFileID)), 
         catch_prop = nByReg / nSampled) %>% 
  ungroup() %>%  
  select(strata, statArea, week, month_week, regName, catch_prop, nByReg, 
         nSampled) %>%
  distinct() %>% 
  mutate(statArea = fct_relevel(as.factor(statArea), "123", after = Inf),
         regName = fct_relevel(regName, "North/Central BC", "WCVI",
                               "Fraser River", "SOG", "Puget Sound", 
                               "Washington Coast", "Oregon/California"))

weekly_reg_plot <- ggplot(weekly_comp) +
  geom_bar(aes(x = as.factor(week), y = catch_prop, fill = regName),
           stat = "identity") +
  scale_fill_viridis_d(option = "B") +
  labs(x = "", y = "Composition", fill = "Region of Origin") +
  facet_wrap(~statArea) + 
  ggsidekick::theme_sleek()

pdf(here::here("figs", "stockComp", "weekly_comp_agg.pdf"), height = 6,
    width = 10)
weekly_reg_plot
dev.off()

## Calculate stock proportion within Fraser assuming perfect assignment for each
# stat area-week and lumping all years; focus only 123/124 for now
# Identify Fraser IDs from above
fraser_ids <- max_prob %>% 
  filter(regName == "Fraser River") %>% 
  pull(flatFileID)
fraser <- fullGSI %>% 
  filter(flatFileID %in% fraser_ids)
max_prob2 <- fraser %>% 
  group_by(flatFileID, Region1Name) %>% 
  mutate(muProb = sum(adjProb)) %>% 
  group_by(flatFileID) %>% 
  mutate(max_assignment = max(muProb),
         strata = paste(statArea, week, sep = "_")) %>% 
  # Remove samples where top stock ID is less than 75% probability
  filter(!muProb < max_assignment, 
         !max_assignment < 0.75,
         statArea %in% c("123", "124"), 
         month %in% c(5, 6, 7, 8)) %>% 
  ungroup() %>% 
  select(-prob, -stock, -adjProb) %>% 
  distinct()

weekly_comp2 <- max_prob2 %>%
  group_by(strata) %>% 
  mutate(nSampled = length(unique(flatFileID))) %>% 
  group_by(strata, Region1Name) %>% 
  mutate(nByMU = length(unique(flatFileID)), 
         catch_prop = nByMU / nSampled) %>% 
  ungroup() %>%  
  select(strata, statArea, week, MU = Region1Name, catch_prop, nByMU, 
         nSampled) %>%
  distinct() %>% 
  mutate(statArea = fct_relevel(as.factor(statArea), "123", after = Inf))

samp_size <- weekly_comp2 %>%
  mutate(week = as.factor(week)) %>% 
  select(statArea, week, nSampled) %>% 
  distinct()

weekly_mu_plot <- ggplot(weekly_comp2) +
  geom_bar(aes(x = as.factor(week), y = catch_prop, fill = MU),
           stat = "identity") +
  scale_fill_viridis_d() +
  labs(x = "Stat Week", y = "Composition", fill = "Management Unit") +
  facet_wrap(~statArea) + 
  ggsidekick::theme_sleek() +
  ylim(0, 1.05) +
  geom_text(data = samp_size, aes(x = week, y= 1.025, 
                                  label = nSampled),
            colour="grey20", size=3)

pdf(here::here("figs", "stockComp", "weekly_comp_FraserMU.pdf"), height = 6,
    width = 10)
weekly_mu_plot
dev.off()
