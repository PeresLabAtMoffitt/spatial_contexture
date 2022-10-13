# Import package
library(targets)
library(tidyverse)

# Load data
tar_load(KNN_ROI_overall_bivariate)
tar_load(KNN_ROI_overall_univariate)
tar_load(C_ROI_overall_bivariate)
tar_load(C_ROI_overall_univariate)

KNN_ROI_overall_bivariate <- KNN_ROI_overall_bivariate %>% 
  janitor::clean_names()
# Data exploratory

table(KNN_ROI_overall_bivariate$r)
# Why do we not have more r? Should have 10, 20, 30, etc
# Why do we see so little cells CD3, it should be more than CD11

is.na(KNN_ROI_overall_bivariate$degree_of_clustering_permutation)

# Choose the radius (r) to use for estimating the spatial measures (i.e. size of neighborhood) 
# should be based on the scale of clustering of interest (i.e. a small value of r will have 
# clustering assess for small neighborhoods while a larger value of r would determine the level 
# of clustering based on large-sized neighbors)
KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation))+
  geom_point()

KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation, color= image_location))+
  geom_line()+ 
  theme(legend.position = "none")

KNN_ROI_overall_bivariate %>% 
  filter(!is.na(degree_of_clustering_permutation)) %>% 
  pivot_wider(id_cols = c(image_location, anchor, counted), 
              names_from = r, 
              values_from = degree_of_clustering_permutation) %>% 
  mutate(count_change = round(`100` - `50`, 2)) %>% 
  ggplot(aes(x= count_change, fill= anchor))+
  geom_bar(stat = "count")+
  ggtitle("difference between radius 100 - 50")+
  theme_classic()+
  facet_wrap(.~ counted)

KNN_ROI_overall_bivariate %>% 
  filter(!is.na(degree_of_clustering_permutation)) %>% 
  group_by(image_location, anchor, counted) %>% 
  mutate(n = n()) %>% 
  filter(n == 2) %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation, color= image_location, group = image_location))+
  geom_line()+
  ggtitle("Ex: CD8 better radius at 50")+
  theme_classic()+
  facet_wrap(counted ~ anchor)+ 
  theme(legend.position = "none")

KNN_ROI_overall_bivariate[23:24,] %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation), color= image_location)+
  geom_line()
KNN_ROI_overall_bivariate[239:240,] %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation), color= image_location)+
  geom_line()

KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= as.factor(r), y= degree_of_clustering_permutation))+
  geom_boxplot()
