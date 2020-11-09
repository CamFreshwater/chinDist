## Principal Tensor Analysis 
# Oct. 6 2020
# Explore using PTA to identify clusters of stocks with similar characteristics
# using predicted abundance from splines

library(tidyverse)
library(ade4)
library(PTAk)

pred_dat <- readRDS(here::here("generated_data",
                               "combined_model_predictions.RDS"))

# Convert predictions ----------------------------------------------------------

# 1. scale estimates of abundance WITHIN each region to account for differences
# in catchability, then scale WITHIN species to generate seasonal anomalies,
# then split into a list by region

scaled_pred_dat <- pred_dat %>% 
  filter(dataset %in% c("gsi_troll_pst", "gsi_sport_pst")) %>% 
  select(dataset, comp_pred_ci) %>% 
  mutate(scaled_pred = map(comp_pred_ci, function (x) {
    x %>% 
      select(month_n, region, stock, comp_abund_est) %>% 
      mutate(scale_abund = scale(comp_abund_est, center = TRUE, 
                                 scale = TRUE)[ , 1]) 
  })) %>% 
  select(-comp_pred_ci) %>% 
  unnest(scaled_pred) %>% 
  select(-dataset, -comp_abund_est) %>%
  # either filter to remove months or filter fewer months and drop some regions
  filter(!month_n < 6,
         !month_n > 9,
         !region == "QCaJS") %>% 
  droplevels()

ggplot(scaled_pred_dat) +
  geom_boxplot(aes(x = stock, y = scale_abund))

pred_list <- scaled_pred_dat %>% 
  #scale within a stock to generate anomalies
  group_by(stock) %>% 
  mutate(scale_abund = scale(scale_abund, center = TRUE, scale = TRUE)[ , 1],
         month = as.character(month_n)) %>%
  select(-month_n) %>% 
  ungroup() %>% 
  pivot_wider(names_from = month,
              values_from = scale_abund) %>% 
  split(., .$region) 

# 2. convert list of dataframes into array
months = colnames(pred_list[[1]])[-c(1,2)]
stocks = unique(scaled_pred_dat$stock)
regions = unique(scaled_pred_dat$region)

pred_array <- array(0, 
                    dim = c(length(stocks), length(months), length(regions)),
                    dimnames = list(stocks, months, regions))
for (i in seq_along(regions)) {
  dum_mat <- pred_list[[i]] %>% 
    select(-stock, -region) %>% 
    as.matrix()
  pred_array[, , i] <- dum_mat
}


# Fit PTA ----------------------------------------------------------------------

pta <- PTA3(pred_array, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta, testvar = 0)

out <- !substr(pta[[3]]$vsnam, 1, 1) == "*"
gct <- (pta[[3]]$pct * pta[[3]]$ssX / pta[[3]]$ssX[1])[out]
barplot(sort(gct, decreasing = TRUE), xlab="PT",
        ylab="Percentage of variance")
#sharp reduction after 3rd PT

# identify tensor product names based on number selected
tp_keep <- which((pta[[3]]$pct * pta[[3]]$ssX / pta[[3]]$ssX[1]) %in% 
                   sort(gct, decreasing = TRUE)[1:3])
pta[[3]]$vsnam[tp_keep]

# plot (mod refers dimensions, nb to the components)
plot(pta, mod=c(1, 2, 3), nb1 = 1, nb2 = 6, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 6, lengthlabels = 6)
plot(pta, mod=c(1, 2, 3), nb1 = 1, nb2 = 7, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 7, lengthlabels = 6)

# Conduct clustering as before
coo <- t(pta[[1]]$v[tp_keep, ])
labkeep <- paste0(pta[[3]]$vsnam[tp_keep], " - ", round((100 * (pta[[3]]$d[tp_keep])^2)/pta[[3]]$ssX[1],1), "%")
rownames(coo) <- stocks

#1. Compute the distance between species
dist1 <- dist(coo, method = "euclidean")

#2. Build a tree with Ward linkage
den <- hclust(dist1,method = "ward.D2")

#3. Plot the dendogram
plot(den, hang=-1, ax = T, ann=F, xlab="", sub="")

#Choose the number of clusters
nclust <- 3
rect.hclust(den, k=nclust, border=rainbow(nclust)[c(3,2,1)])
