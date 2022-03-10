# Import library
library(tidyverse)

################################################# I ### Load data
path <- fs::path("","Volumes","Peres_Research", "K99_R00", "Image analysis data", "Spatial Data")

ROI_tumor <- 
  read_csv(
    paste0(path,
           "/Nearest Neighbor/ROI_tumor_bivariate.csv"))

list.files(path = paste0(path, "/Nearest Neighbor/"),
           pattern = "*.csv", 
           full.names = T) %>% 
  map_df(~read_csv(., col_types = cols(.default = "c"))) 
