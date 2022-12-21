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
  ungroup() %>% 
  mutate(image_tag = str_match(image_location, "Analysis Images.(.*?).tif")[,2],
         image_tag = str_replace(image_tag, "16-", "16")) %>% 
  mutate(suid = str_match(image_tag,
                          "(L.Peres_P1_OV|L.Peres_P1_)([:digit:]*)")[,3],
         .before = 1) %>%
  select(-image_location)

KNN_ROI_overall_bivariate <- KNN_ROI_overall_bivariate %>% 
  filter(r == 50 & anchor != "DAPI (DAPI) Positive") %>% 
  select(-r)


# ICC
df <- KNN_ROI_overall_bivariate %>% 
  as.data.frame(.)

"theoretical_csr"                  "permuted_g"                      
[7] "observed_g"                       "degree_of_clustering_permutation" "degree_of_clustering_theoretical"

library(psych)

fct_icc <- function(data) {
  ICC_data <- data.frame(matrix(nrow = 1, ncol = 0))
  lb_data <- data.frame(matrix(nrow = 1, ncol = 0))
  up_data <- data.frame(matrix(nrow = 1, ncol = 0))
  for (i in 1:length(colnames(data))) {
    rad <- data %>% select(suid, anchor, counted)
    
    if (class(data[, i]) == "numeric" |
        class(data[, i]) == "integer") {
      ICC_df <- cbind(rad, value = data[, i])
      ICC_df <- ICC_df %>%
        mutate(Slide = "Slide0") %>%
        group_by(suid) %>%
        mutate(n = row_number(suid)) %>%
        ungroup() %>%
        unite(slide_id, Slide:n, sep = "", remove = TRUE, na.rm = TRUE) %>%
        pivot_wider(id_cols = -c(anchor, counted), 
                    names_from = slide_id, values_from = value) %>%
        select(c(starts_with("slide")))
      
      ICC <- ICC(ICC_df)$results[4, 2]
      ICC_data <- cbind(ICC_data, ICC)
      ICC <- ICC(ICC_df)$results[4, 7]
      lb_data <- cbind(lb_data, ICC)
      ICC <- ICC(ICC_df)$results[4, 8]
      up_data <- cbind(up_data, ICC)
    }
  }
  ICCC_data <- bind_rows(ICC_data, lb_data, up_data)
  colnames(ICCC_data) <- colnames(data)[6:ncol(data)-1]
  ICCC_data <- as.data.frame(t(ICCC_data)) %>%
    mutate(ICC_lb_up = paste(round(V1, 2), " (", round(V2, 2), ", ", round(V3, 2), ")", sep = "")) %>%
    select(ICC_lb_up)
}

dat <- fct_icc(df)
dat









########### NEXT can summarize########### 
# Let's choose 50 for now

KNN_ROI_overall_bivariate <- KNN_ROI_overall_bivariate %>% 
  group_by(suid, anchor, counted) %>% 
  summarize(theoretical_csr = mean(theoretical_csr, na.rm=TRUE),
            permuted_g = mean(permuted_g, na.rm=TRUE),
            observed_g = mean(observed_g, na.rm=TRUE),
            degree_of_clustering_permutation = mean(degree_of_clustering_permutation, na.rm=TRUE),
            degree_of_clustering_theoretical = mean(degree_of_clustering_theoretical, na.rm=TRUE)) %>% 
  ungroup()




# Data exploratory

table(KNN_ROI_overall_bivariate$r)
# Why do we not have more r? Should have 10, 20, 30, etc
# Why do we see so little cells CD3, it should be more than CD11

is.na(KNN_ROI_overall_bivariate$degree_of_clustering_permutation)

########### CHOOSE R ########### 
# Choose the radius (r) to use for estimating the spatial measures (i.e. size of neighborhood) 
# should be based on the scale of clustering of interest (i.e. a small value of r will have 
# clustering assess for small neighborhoods while a larger value of r would determine the level 
# of clustering based on large-sized neighbors)
################################
KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation))+
  geom_point()

KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation, color= image_tag))+
  geom_line()+ 
  theme(legend.position = "none")

KNN_ROI_overall_bivariate %>% 
  filter(!is.na(degree_of_clustering_permutation)) %>% 
  pivot_wider(id_cols = c(image_tag, anchor, counted), 
              names_from = r, 
              values_from = degree_of_clustering_permutation) %>% 
  mutate(count_change = round(abs(`100`) - abs(`50`), 2)) %>% 
  ggplot(aes(x= count_change, fill= anchor))+
  geom_bar(stat = "count")+
  ggtitle("difference between radius 100 - 50, if more cells at 100 -> postive value")+ # Why do we have more at 50
  theme_classic()+
  facet_wrap(.~ counted)

KNN_ROI_overall_bivariate %>% 
  filter(!is.na(degree_of_clustering_permutation)) %>% 
  group_by(image_tag, anchor, counted) %>% 
  mutate(n = n()) %>% 
  filter(n == 2) %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation, color= image_tag, group = image_tag))+
  geom_line()+
  ggtitle("Ex: CD8 better radius at 50")+
  theme_classic()+
  facet_wrap(counted ~ anchor)+ 
  theme(legend.position = "none")

KNN_ROI_overall_bivariate[23:24,] %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation), color= image_tag)+
  geom_line()
KNN_ROI_overall_bivariate[239:240,] %>% 
  ggplot(aes(x= r, y= degree_of_clustering_permutation), color= image_tag)+
  geom_line()

KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= as.factor(r), y= degree_of_clustering_permutation))+
  geom_boxplot()








KNN_ROI_overall_bivariate %>% 
  ggplot(aes(x= counted, y= degree_of_clustering_permutation, color= counted))+
  geom_boxplot()+
  ggtitle("Ex: CD8 better radius at 50")+
  theme_classic()+
  facet_wrap(. ~ anchor)+
  coord_flip()+ 
  theme(legend.position = "none")




################################








