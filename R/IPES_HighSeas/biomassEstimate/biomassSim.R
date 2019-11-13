### Simulate biomass estimates to quantify relative benefits of increasing 
# sampling effort within a strate

###### IMPORT DATA AND CLEAN
library(tidyverse)
library(ggplot2)
fish <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                            "Biomass_Estimates.csv"), stringsAsFactors = F) 

volswept <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                                "VolumeSwept_ByStratum_DayNight.csv"), 
                     stringsAsFactors = FALSE)
haul <- read.csv(here::here("data", "highSeas", "ipesBiomass", 
                            "Trawl_Tows.csv"), stringsAsFactors = F) %>%  
  merge(., volswept,by = c("TRIP_YEAR","STRATUM","DayNight"), 
        all.x = TRUE) %>% 
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER))


df <- fish %>% 
  dplyr::rename(SPECIES_CODE_orig = SPECIES_CODE) %>%
  mutate(SPECIES_CODE = as.integer(str_sub(SPECIES_CODE_orig, 1, 3))) %>%
  mutate(yearTow = paste0(TRIP_YEAR, "-", TOW_NUMBER),
         CatchWt = 
           if_else(SPECIES_CODE == 96, 
                   TotalCatchWtCPUE_kg_km3, JuvCatchWtCPUE_kg_km3) / 1000) %>% 
  select(TRIP_YEAR, STRATUM, yearTow, SPECIES_CODE, CatchWt)

# add zeros to catch data and merge with hauls
cpue_df <- haul %>% 
  select(yearTow, TRIP_YEAR, STRATUM, DayNight) %>% 
  left_join(., df %>% select(yearTow, SPECIES_CODE), by = "yearTow") %>% 
  expand(., nesting(yearTow, TRIP_YEAR, STRATUM, DayNight), SPECIES_CODE) %>% 
  filter(!is.na(SPECIES_CODE)) %>% 
  left_join(., df, by = c("yearTow", "TRIP_YEAR", "STRATUM", "SPECIES_CODE")) %>% 
  complete(yearTow, SPECIES_CODE, fill = list(CatchWt = 0)) %>% 
  mutate(nonZero = case_when(
    CatchWt > 0 ~ 1,
    CatchWt == 0 ~ 0
  )) 


#1. Use gamma hurdle models to estimate parameters describing data then simulate
modfits <- cpue_df %>% 
  group_by(SPECIES_CODE) %>% 
  nest()  %>% 
  #don't account for year/strata now, but could as fixed or mixed effects
  mutate(m1 = map(data, ~glm(nonZero ~ STRATUM, data = ., 
                             family = binomial(link = logit))),
         m2 = map(data, 
                  ~glm(CatchWt ~ STRATUM, data = subset(., nonZero == 1), 
                       family = Gamma(link = log)))
         ) %>% 
  mutate(tidy1 = map(m1, broom::tidy),
         tidy2 = map(m2, broom::tidy),
         summ = map(m2, summary),
         gamDisp = map_dbl(summ, "dispersion"))

pvals <- modfits %>% 
  unnest(tidy1, tidy2) %>% 
  select(species = SPECIES_CODE, term, binMu = estimate,  
         binP = p.value, gamMu = estimate1, gamP = p.value1) %>% 
  pivot_longer(-c(species, term), names_to = "outputValue", 
               values_to = "estimate") %>%
  mutate(
    model = case_when(
      grepl("gam", outputValue) ~ "gamma",
      TRUE ~ "binomial"),
    value = case_when(
      grepl("Mu", outputValue) ~ "mu",
      grepl("P", outputValue) ~ "p-value"
    )) %>%
  select(-outputValue)
  
pvals %>% 
  filter(term == "STRATUM", 
         value == "p-value")
#generally strata effects are significant for the binomial model, but less 
#consistent for gamma

modOut <- modfits
  unnest(tidy1, tidy2, gamDisp) %>%
  select(species = SPECIES_CODE, data, binMu = estimate, binMuSig = std.error, 
         gamMu = estimate1, gamMuSig = std.error1, gamDisp) %>% 
  pivot_longer(-c(species, data), names_to = "parameter", 
               values_to = "estimate") %>% 
  mutate(
    model = case_when(
      grepl("gam", parameter) ~ "gamma",
      TRUE ~ "binomial"),
    parameter = case_when(
      grepl("MuSig", parameter) ~ "sigmaMu",
      grepl("Mu", parameter) ~ "mu",
      grepl("Disp", parameter) ~ "dispersion"
    )) %>% 
  arrange(species, model, parameter)


  
  
# Best practice would be to simulate data by specifying trials; each trial 
# generate an estimate of phi and beta (rate) by drawing from normal 
# distribution describing mu and se of parameter estimates. Then use these values
# to generate 1000 fishing events (i.e. 0/1 and estimate of catch if non-0). 
# Eventually this should account for strata and year effects as well.

# simpler example with one species  
dum <- cpue_df %>% 
  filter(SPECIES_CODE == "124")
m1 <- glm(nonZero ~ STRATUM, data = dum, family = binomial(link = logit))
m2 <- glm(CatchWt ~ STRATUM, data = subset(dum, nonZero == 1), 
          family = Gamma(link = log))

newDat <- data.frame(STRATUM = unique(dum$STRATUM))
#binomial predictions
phiStrata <- predict(m1, newdata=newDat, type='response') # on scale of response
#gamma predictions
muStrata <- predict(m2, newdata = newDat, type = 'response')
#gamma distribution parameters
gammaCoefs <- data.frame(strata = unique(dum$STRATUM),
                         shape = 1 / summary(m2)$dispersion) %>%
  mutate(rate = shape / muStrata)

# generate random draws for the binomial catch
M <- 10 #number of trials
N <- 10 #number of draws per trial
trialList  <- data.frame(strata = rep(unique(dum$STRATUM), each = M),
                   trial = rep(seq(1, M, by = 1), 
                              times = length(unique(dum$STRATUM))),
                   #binomial coefficients first
                   phi = rep(phiStrata, each = M)) %>% 
  left_join(., gammaCoefs, by = "strata") %>% 
  split(., .$trial) %>% 
  lapply(., function(x) {
    #simulateTrials func here
    draw <- seq(1, N, by = 1)
    x %>% 
      expand()
  })
  
simulateTrials <- function(dattt, N) {
  dattt <- trialList[[1]] 
  draw <- seq(1, N, by = 1)
  # Sample binomial distribution within strata/trial N draw times
  tempBinomial <- dattt %>% 
    expand(nesting(strata, trial, phi, shape, rate), draw) %>% 
    group_by(strata) %>% 
    #don't account for year/strata now, but could as fixed or mixed effects
    mutate(nonZero = rbinom(N, size = 1, prob = phi)) 
  
  # Sample gamma distribution within strata/trial for nonzero trials within 
  # a strata
  tempGamma <- tempBinomial %>% 
    filter(nonZero == "1") %>% 
    group_by(strata) %>%
    mutate(towsWCatch = length(nonZero)) %>% 
    mutate(catchWt = rgamma(n =  towsWCatch, shape = shape, rate = rate)) %>% 
    select(strata, draw, catchWt)
  catchOut <- tempBinomial %>% 
    left_join(., tempGamma, by = c("strata", "draw")) %>% 
    replace_na(., list(catchWt = 0))
}
 



catchOut <- dum2 %>% 
  left_join(., gammaDum, by = c("strata", "draw")) %>% 
  replace_na(., list(catchWt = 0))


ggplot(catchOut %>% filter(nonZero == "1")) +
  geom_histogram(aes(x = catchWt)) +
  facet_wrap(~strata)

ggplot(cpue_df %>% filter(SPECIES_CODE == "124",
                          nonZero == "1")) +
  geom_histogram(aes(x = CatchWt)) +
  facet_wrap(~STRATUM)
  
#2. Sample from distribution at different levels
#3. Calculate CV following biomass extrapolotation


set.seed(999)
N <- 100
x <- runif(N, -1, 1)
a <- 0.5
b <- 1.2
y_true <- exp(a + b *x)
shape <- 10
y <- rgamma(N, rate = shape / y_true, shape = shape)
m_glm <- glm(y ~ x, family = Gamma(link = "log"))
summary(m_glm)

(exp(0.47255)^2) / (exp(0.02925)^2)
1 / 0.08556405

summary(m_glm)$dispersion/coef(m_glm)[1]
#shape parameter
alpha <- 1 / summary(m_glm)$dispersion
#rate parameter (shape / mean)
beta <- summary(m_glm)$dispersion / coef(m_glm)[1] 


x <- rnorm(100) 
y <- rpois(rep(1,100), exp(x)) ## poisson regression with slope=1
## fit model
m1 <- glm(y ~ x,family=poisson)
new.data <- data.frame(x=seq(-3,3,.1))
mu.y <- predict(m1, newdata=new.data, type='response')
sim.y <- replicate(2, rpois(rep(1, length(mu.y)), mu.y))
