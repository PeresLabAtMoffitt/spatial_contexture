# Import package
library(targets)
library(tidyverse)

# Load data
tar_load(KNN_ROI_overall_bivariate)
tar_load(KNN_ROI_overall_univariate)
tar_load(C_ROI_overall_bivariate)
tar_load(C_ROI_overall_univariate)

# Data exploratory

table(KNN_ROI_overall_bivariate$r)
# Why do we not have more r? Should have 10, 20, 30, etc

is.na(KNN_ROI_overall_bivariate$`Degree of Clustering Theoretical`)

# Choose the radius (r) to use for estimating the spatial measures (i.e. size of neighborhood) 
# should be based on the scale of clustering of interest (i.e. a small value of r will have 
# clustering assess for small neighborhoods while a larger value of r would determine the level 
# of clustering based on large-sized neighbors)
KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= r, y= `Degree of Clustering Theoretical`))+
  geom_point()

KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= r, y= `Degree of Clustering Theoretical`), color= `Image Location`)+
  geom_line()

KNN_ROI_overall_bivariate[23:24,] %>% 
  ggplot(aes(x= r, y= `Degree of Clustering Theoretical`), color= `Image Location`)+
  geom_line()
KNN_ROI_overall_bivariate[239:240,] %>% 
  ggplot(aes(x= r, y= `Degree of Clustering Theoretical`), color= `Image Location`)+
  geom_line()

KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= as.factor(r), y= `Degree of Clustering Theoretical`))+
  geom_boxplot()
