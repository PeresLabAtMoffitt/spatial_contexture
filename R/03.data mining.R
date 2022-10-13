# Import package
library(targets)
library(tidyverse)

# Load data
tar_load(KNN_ROI_overall_bivariate)
tar_load(KNN_ROI_overall_univariate)
tar_load(C_ROI_overall_bivariate)
tar_load(C_ROI_overall_univariate)

KNN_ROI_overall_bivariate <- KNN_ROI_overall_bivariate %>% 
  janitor::clean_names() %>% 
  filter(r == 50 & 
           !str_detect(counted, "DAPI") & 
           !str_detect(anchor, "DAPI"))  ##### But what if I want to see how many CD3 arounf 1 tumor cell?

# Data mining
KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= counted, y= permuted_g))+
  geom_boxplot()+
  ggtitle("G for counted cell type around each anchor cell type")+
  theme_classic()+
  coord_flip()+
  facet_wrap(.~ anchor)

KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= counted, y= degree_of_clustering_permutation))+
  geom_boxplot()+
  ggtitle("clustering for counted cell type around each anchor cell type")+
  theme_classic()+
  coord_flip()+
  facet_wrap(.~ anchor)
