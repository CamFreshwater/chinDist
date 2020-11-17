# Estimating BC Chinook salmon spatio-temporal distributions using fisheries and genetic stock composition data

This repository contains code associated with the paper:

Freshwater, C., S.C. Anderon, T.D. Beacham, W. Luedke, C. Wor, and J. King. 2020. An integrated model for seasonal changes in stock composition and abundance with an application to Chinook salmon in British Columbia. 

## Dependencies

In order to run the included code the following R packages are required:

``` r
pkgs <- c("tidyverse", "ggplot2", "TMB", "ggplot2", "grid", "gridExtra", "mgcv",
"scales", "here")
install.packages(pkgs)
devtools::install_github("seananderson/ggsidekick")
devtools::install_github("pbs-assess/stockseasonr")
```
