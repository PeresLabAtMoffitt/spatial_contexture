# Import library
library(tidyverse)


# tar_load(K_ROI_tumor_bivariate)
################################################# II ### Data cleaning

# script for K data. will create function later


bivariate_k <- K_ROI_tumor_bivariate %>% 
  # Bind compartment data
  mutate(compartment = "Tumor") %>% 
  bind_rows(., K_ROI_overall_bivariate %>%
              mutate(compartment = "Total")) %>%
  bind_rows(., K_ROI_stroma_bivariate %>%
              mutate(compartment = "Stroma")) %>%
  # filter the r value we will use on our analysis
  filter(r %in% c(20:30)) %>% 
  # create average of K value
  group_by(suid, compartment, anchor, counted, r) %>% ##################### Add periph, intra region-------
  summarize(average_CSR = mean(`Theoretical CSR`, na.rm = TRUE),
            average_perm_k = mean(`Permuted K`, na.rm = TRUE),
            average_obs_k = mean(`Observed K`, na.rm = TRUE))


write_rds(bivariate_k, "bivariate_k.rds")
