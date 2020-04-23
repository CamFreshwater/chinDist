# Function to clean composition data (either from GSI or CWT) including 
# infilling zero values and pivoting to wide format
# Should be in long format (e.g. reg3RollUpCatchProb.RDS)

clean_comp <- function(comp, month_range = c(1, 12), check_tables = FALSE) {
  comp_trim <- comp %>% 
    filter(!month_n < month_range[1],
           !month_n > month_range[2]) %>% 
    droplevels() %>% 
    dplyr::select(id, statArea, year, month, month_n, season, 
                  regName, pres, catchReg) 
  
  # multinomial model predictions are a little wonky w/ missing values - this is 
  # one infilling option
  # dummy dataset to replace missing values 
  stock_names <- unique(comp_trim$regName)
  #retain only values from the dummy set that are not already present
  missing_sample_events <- expand.grid(
    id = NA,
    statArea = NA,
    year = NA,
    month = unique(comp_trim$month),
    month_n = NA,
    season = NA,
    regName = stock_names,
    pres = 1,
    catchReg = NA
  ) %>%
    anti_join(., comp_trim,
              by = c("month", "regName", "pres")
    ) %>%
    select(-regName) %>%
    distinct()
  #add random subset of yrs to avoid overparameterizing model
  rand_yrs <- sample(unique(comp_trim$year), size = nrow(missing_sample_events),
                     replace = TRUE)
  rand_catch_reg <- sample(unique(comp_trim$catchReg), 
                           size = nrow(missing_sample_events),
                           replace = TRUE)
  missing_sample_events$year <- rand_yrs
  missing_sample_events$catchReg <- rand_catch_reg
  
  #duplicate for all stocks (to balance the addition)
  infill_data <- do.call("rbind", replicate(length(stock_names),
                                            missing_sample_events,
                                            simplify = FALSE)) %>%
    mutate(regName = rep(stock_names, each = nrow(missing_sample_events))) %>%
    select(id:season, regName, pres, catchReg)
  
  # combine and check for no zeros
  temp <- comp_trim %>%
    rbind(., infill_data)
  tab_in <- table(comp_trim$regName, comp_trim$month, comp_trim$catchReg)
  tab_out <- table(temp$regName, temp$month, temp$catchReg)
  
  # spread
  comp_out <- temp %>% 
    mutate(dummy_id = seq(from = 1, to = nrow(.), by = 1)) %>% 
    pivot_wider(., names_from = regName, values_from = pres) %>%
    mutate_if(is.numeric, ~replace_na(., 0)) %>%
    droplevels() 
  
  # comp_out %>% filter(month == "10", `CR-sp&su` > 0) %>% select(`CR-sp&su`)
  # temp %>% filter(month == "1", regName == "Other") %>% 
  #   select(regName, month, pres)
  # 
  ifelse(check_tables == FALSE, 
         return(comp_out), 
         return(list(data = comp_out, long_data = temp,
                     tables = list(original = tab_in, infilled = tab_out)))
  )
}