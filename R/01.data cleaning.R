# Import package
library(targets)
library(tidyverse)

# Load data
tar_load(KNN_ROI1_overall_univariate)
# tar_load(C_ROI1_overall_univariate)
tar_load(C_ROI1_overall_univariate)
# tar_load(C_ROI1_overall_univariate)


# Clean names and missing
# + Extract Ids
clean_data <- function(data) {
  data <- data %>% 
    janitor::clean_names() %>% 
    ungroup() %>%
    mutate(across(where(is.numeric), ~ na_if(., NaN))) %>% 
    mutate(image_tag = str_match(image_location, "Analysis Images.(.*?)$")[,2],
           image_tag = str_replace(image_tag, "16-", "16")) %>% 
    mutate(suid = str_match(image_tag,
                            "(Peres_P1_OV|Peres_P1_)([:digit:]*)")[,3],
           .before = 1) %>%
    select(-image_location)
}

KNN_ROI1_overall_univariate <- clean_data(KNN_ROI1_overall_univariate)
write_rds(KNN_ROI1_overall_univariate, "knn_ROI1_overall_univariate.rds")

C_ROI1_overall_univariate <- clean_data(C_ROI1_overall_univariate)
write_rds(C_ROI1_overall_univariate, "c_ROI1_overall_univariate.rds")








