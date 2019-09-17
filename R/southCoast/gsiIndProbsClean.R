## Clean individual genetic data
# Data from west coast Vancouver Island commercial troll fishery
# Sep 14, 2019

library(tidyverse)
library(ggplot2)

datRaw <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbs_trim.csv"), 
                stringsAsFactors = FALSE)

# pull ID column, split and rename
id <- data.frame(t(matrix(unlist(strsplit(as.vector(datRaw$sampleID), 
                                          split = " ")),
                        nrow = 5, ncol = length(datRaw$sampleID)))) %>% 
  rename("statArea" = X1, "year" = X2, "gear" = X3, "jDay" = X4, 
         "fishNum" = X5) %>% 
  mutate(statArea = case_when(
            statArea %in% c("Area023", "Area23", "Area_23") ~ "23",
            statArea %in% c("Area123", "Area123SWVI", "Area123Comm") ~ "123",
            statArea %in% c("Area124", "Area124SWVI", "Area124Comm") ~ "124",
            statArea %in% c("Area125", "Area125NWVI") ~ "125",
            statArea %in% c("Area126", "Area126NWVI") ~ "126",
            statArea %in% c("Area127", "Area127NWVI") ~ "127",
            statArea %in% c("Area026") ~ "26",
            statArea %in% c("Area24", "Area_24", "Area124_24", 
                            "Area_24xgill") ~ "24",
            TRUE ~ as.character(statArea)),
         abbYear = sapply(strsplit(as.character(year), '[()]'), 
                       function(x) (x)[2]),
         year = paste("20", abbYear, sep = ""),
         date = as.Date(as.numeric(as.character(jDay)), 
                        origin = as.Date(paste(year, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")))

dat <- cbind(id, datRaw) %>% 
  select(-sampleID, -abbYear) %>% 
  # for now remove ~150 fish that can't be assigned to an individual stat area
  filter(!statArea %in% c("Area123-124", "Area125-126", "Area126-127", 
                          "Area124_24"))

temp <- dat %>% 
  group_by(statArea, year, month) %>% 
  tally()

summDat <- expand.grid(unique(dat$statArea), unique(dat$year), unique(dat$month)) %>%
  rename("statArea" = Var1, "year" = Var2, "month" = Var3) %>% 
  left_join(., temp, by = c("statArea", "year", "month")) %>% 
  replace_na(list(n = 0))

