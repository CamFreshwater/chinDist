## Clean CWT data
# Import and clean CWT data from MRP extractor:
# http://dfbcv8lwvast012.ent.dfo-mpo.ca/DataExtractor/#/Query
# then generate equivalent data to GSI input data to pass to multinomial or 
# combined tmb models.
# Note that previous version used RMIS data, but could not be easily assigned
# to PFMAs, which is necessary for assigning to current catch regions.
# MRP data is equivalent, but will lack non-DFO recoveryies.#
# Change made June 22, 2020

library(tidyverse)

## ADD STOCKS TO RECOVERY DATA -------------------------------------------------

#all fisheries recoveries from mrp extractor 
#technically this dataset could be used for the analysis, but would require
#additional cleaning for now, use to generate key for meaningful PFMA names 
# and the RMIS recoveries (above)
tag_recov <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                                 "mrp_recoveries.csv"), 
                      stringsAsFactors = FALSE) %>% 
  #remove samples without tags
  filter(!Tagcode == "No Tag") %>% 
  mutate(nchar_tag = nchar(Tagcode)) %>% 
  #remove samples with non-standard tags 
  filter(!nchar_tag > 6,
         !Tagcode %in% c("53397", "54380", "181986")) %>% 
  mutate(Tagcode = stringr::str_pad(Tagcode, 6, pad = "0"))
  
# export tag codes to query MRP extractor (or rmis)
tag_codes <- unique(tag_recov$Tagcode)
# write.table(tag_codes, here::here("data", "gsiCatchData", "commTroll",
#                                   "cwt_tag_code_mrp.txt"),
#             row.names = FALSE, col.names = FALSE)

# import generated tag_code key (from rmis standard reporting query)
mrp_key_raw <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                                    "cwt_tag_key_mrp.txt"),
                         stringsAsFactors = FALSE)

# generate trimmed clean tag key to pass to stockKey repo for merging
# clean_key <- mrp_key_raw %>%
#   mutate(
#     #when stock name missing replace with hatchery name
#     stock_location_name = case_when(
#       is.na(stock_location_name) ~ hatchery_location_name,
#       TRUE ~ stock_location_name
#     )
#   ) %>%
#   select(release = release_location_name,
#          stock = stock_location_name,
#          rmis_region = release_location_rmis_region,
#          basin = release_location_rmis_basin,
#          state = release_location_state) %>%
#   distinct()
# saveRDS(clean_key, here::here("data", "stockKeys", "cwt_stock_key_out.RDS"))

# import completed stock_list
stock_key <- readRDS(here::here("data", "stockKeys",
                             "finalStockList_June2020.rds")) 

# add regional aggregates to key
key_final <- mrp_key_raw %>%
  select(tag_code = tag_code_or_release_id,
         stock = stock_location_name, 
         basin = release_location_rmis_basin) %>%
  mutate(tag_code = stringr::str_pad(tag_code, 6, pad = "0"),
         stock = toupper(stock),
         stock = case_when(
           basin == "CLEA" ~ "SNAKE R FALL",
           TRUE ~ stock
         )) %>% 
  select(-basin) %>% 
  left_join(., stock_key,
            by = c("stock"
                   )) %>% 
  filter(!is.na(Region1Name))

#merge stock ID and recovery location keys to cwt recoveries
tag_recov1 <- tag_recov %>% 
  select(tag_code = Tagcode, recovery_agency = Recovery.Agency.Code,
         year = Recovery.Year, month = Recovery.Month, 
         Catch.Region.Name:MRP.Area.Name, observed = Observed.Number, 
         est = Estimated.Number, date = Recovery.Date) %>% 
  left_join(., key_final, by = "tag_code") %>% 
  # convoluted extraction process to identify PFMAs then roll up to catch regions
  mutate(
    area_code_trim = case_when(
      grepl("A", MRP.Area.Code) | grepl("B", MRP.Area.Code) | 
        grepl("M", MRP.Area.Code) ~ 
        stringr::str_sub(MRP.Area.Code, start = 2, end = 3),
      TRUE ~
        stringr::str_sub(MRP.Area.Code, start = 2, end = 4)),
    pfma = case_when(
      area_code_trim == "SPS" ~ "sSoG_composite",
      area_code_trim == "40" ~ "SWVI_composite",
      area_code_trim == "05" ~ "NWVI_composite",
      area_code_trim == "03" ~ "nSoG_composite",
      TRUE ~ stringr::str_remove(area_code_trim, "^0+")
    ),
    pfma_n = as.numeric(pfma),
    region = case_when(
      pfma_n > 124 ~ "NWVI",
      pfma_n < 28 & pfma_n > 24 ~ "NWVI",
      pfma_n %in% c("20", "121", "21") ~ "Juan de Fuca Strait",
      is.na(pfma_n) ~ "Juan de Fuca Strait",
      pfma_n < 125 & pfma_n > 120 ~ "SWVI",
      pfma_n < 25 & pfma_n > 20 ~ "SWVI",
      pfma_n %in% c("14", "15", "16") ~ "N. Strait of Georgia",
      pfma_n %in% c("17", "18", "19", "28", "29") ~ "S. Strait of Georgia",
      pfma_n %in% c("10", "11", "111") ~ "Queen Charlotte Sound",
      pfma_n %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits",
      grepl("SoG", pfma) ~ "Georgia Strait",
      grepl("SWVI", pfma) ~ "SWVI",
      grepl("NWVI", pfma) ~ "NWVI",
      TRUE ~ NA_character_
    ),
    gear = case_when(
      grepl("Sport", Catch.Region.Name) ~ "sport",
      grepl("Troll", Catch.Region.Name) ~ "troll",
      grepl("Taaq", Catch.Region.Name) ~ "troll"
    ) %>% 
      as.factor()
    )  


## FORMAT CWT DATA TO MATCH GSI ------------------------------------------------

rec_gsi <- readRDS(here::here("data", "gsiCatchData", "rec", 
                          "recIndProbsLong.rds"))
comm_gsi <- readRDS(here::here("data", "gsiCatchData", "commTroll", 
                           "wcviIndProbsLong.rds"))

tag_recov_out <- tag_recov1 %>%
  # drop years and regions from cwt data that are missing in gsi
  mutate(
    include = case_when(
      gear == "troll" & region %in% comm_gsi$region & 
        year %in% comm_gsi$year ~ T,
      gear == "sport" & region %in% rec_gsi$region & 
        year %in% rec_gsi$year ~ T,
      TRUE ~ FALSE
    )
  ) %>% 
  filter(include == T,
         #remove small number of tag codes that didn't show up in RMIS
         !tag_code %in% c("0BLANK", "054380", "053397", "184376")) %>%
  # add in relevant helper variables and rename
  mutate(temp_strata = paste(month, region, sep = "_"),
         season_c = case_when(
           month %in% c("12", "1", "2") ~ "w",
           month %in% c("3", "4", "5") ~ "sp",
           month %in% c("6", "7", "8") ~ "su",
           month %in% c("9", "10", "11") ~ "f"
         ),
         season = fct_relevel(season_c, "sp", "su", "f", "w"),
         month_n = as.numeric(month),
         month = as.factor(month_n),
         pres = 1) %>%
  group_by(gear) %>% 
  nest() %>% 
  mutate(subset_data = map(data, function (x) {
    x %>% 
      mutate(year = as.factor(year)) %>% 
      select(id = tag_code, temp_strata, region, area = pfma, year, month,
             date, #gear, 
             pres, season, month_n, area_n = pfma_n, 
             stock:pst_agg)
  })) 

# tt <- tag_recov_out %>% 
#   unnest(., cols = c(subset_data)) %>% 
#   filter(is.na(pst_agg)) %>% 
#   pull(id) %>% 
#   unique()

saveRDS(tag_recov_out, here::here("data", "gsiCatchData", "commTroll",
                             "cwt_recovery_clean.rds"))


tag_recov_out %>% 
  select(-data) %>% 
  unnest(cols = c(subset_data)) %>% 
  group_by(region, pst_agg) %>% 
  tally() %>% 
  print(n = Inf)
