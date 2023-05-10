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
           # keep in case but looks like the mistake was not introduce in the spatial file
           image_tag = str_replace(image_tag, "16-", "16")) %>% 
    filter(!str_detect(image_tag, "Ctrl")) %>% 
    # name changed from before
    mutate(suid = str_match(image_tag,
                            "(Peres_P1_OV|Peres_P1_AACEES |Peres_P1_AACES |Peres_P1_|Peres_P3_)([:digit:]*)")[,3],
           .before = 1) %>%
    # Fix Ids - few Ids are missing a 1 
    mutate(count_number_in_ids = sapply(str_split(suid, ""), length)) %>% 
    mutate(suid = case_when(
      count_number_in_ids == 6         ~ suid,
      count_number_in_ids == 5 &
        str_detect(suid, "4....")      ~ suid,
      count_number_in_ids == 5 &
        str_detect(suid, "3....")      ~ paste0("1", suid)
    )) %>% 
    # Fix image tag
    mutate(image_tag = case_when(
      str_detect(image_tag, "Peres_P1_3") &
        str_detect(suid, "^1")             ~ str_replace(image_tag, "Peres_P1_3", "Peres_P1_13"),
      TRUE   ~ image_tag
    )) %>% 
    select(-c(image_location, count_number_in_ids))
}


KNN_ROI1_overall_univariate <- clean_data(KNN_ROI1_overall_univariate)
write_rds(KNN_ROI1_overall_univariate, "knn_ROI1_overall_univariate.rds")

C_ROI1_overall_univariate <- clean_data(C_ROI1_overall_univariate)
write_rds(C_ROI1_overall_univariate, "c_ROI1_overall_univariate.rds")


# END Cleaning
