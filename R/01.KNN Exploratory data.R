# Import package
library(targets)
library(tidyverse)
library(gtsummary)
library(survival)
library(survminer)

# Load data
tar_load(KNN_ROI_overall_bivariate)
# tar_load(KNN_ROI_overall_univariate)
# tar_load(C_ROI_overall_bivariate)
# tar_load(C_ROI_overall_univariate)

KNN_ROI_overall_bivariate1 <- KNN_ROI_overall_bivariate %>% 
  janitor::clean_names() %>% 
  ungroup() %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN))) #%>% 
  
  
  
KNN_ROI_overall_bivariate2 <- KNN_ROI_overall_bivariate1 %>% 
  mutate(image_tag = str_match(image_location, "Analysis Images.(.*?)$")[,2],
         image_tag = str_replace(image_tag, "16-", "16")) %>% 
  mutate(suid = str_match(image_tag,
                          "(L.Peres_P1_OV|L.Peres_P1_)([:digit:]*)")[,3],
         .before = 1) %>%
  select(-image_location)

markers <- readRDS("~/Documents/GitHub/Peres/spatial_contexture/complete_AACES_NCOCS_batch1_2_07072022.rds")

markers <- markers %>% 
  select(image_tag, suid, annotation, slide_type, site, data_version)

KNN_ROI_overall_bivariate3 <- KNN_ROI_overall_bivariate2 %>% 
  left_join(., markers %>% 
              mutate(is_in_abundance = "Abundance"), 
            by= c("image_tag", "suid"))

# NEED TO FIX IDs in batch2

########################################################################### I ### AUC based on degree of clustering
# area under the Degree of Clustering curve (computed at various r) for choosing the best radius (r)
table(KNN_ROI_overall_bivariate$r)

########### CHOOSE R ########### 
# Choose the radius (r) to use for estimating the spatial measures (i.e. size of neighborhood) 
# should be based on the scale of clustering of interest (i.e. a small value of r will have 
# clustering assess for small neighborhoods while a larger value of r would determine the level 
# of clustering based on large-sized neighbors)
################################
# KNN_ROI_overall_bivariate4 %>% 
#   ggplot(aes(x= r, y= degree_of_clustering_permutation))+
#   geom_point()
# 
# KNN_ROI_overall_bivariate %>% 
#   ggplot(aes(x= r, y= degree_of_clustering_permutation, color= suid))+
#   geom_line()+ 
#   theme(legend.position = "none")
# 
# KNN_ROI_overall_bivariate %>% 
#   filter(!is.na(degree_of_clustering_permutation)) %>% 
#   pivot_wider(id_cols = c(suid, anchor, counted), 
#               names_from = r, 
#               values_from = degree_of_clustering_permutation) %>% 
#   mutate(count_change = round(abs(`100`) - abs(`50`), 2)) %>% 
#   ggplot(aes(x= count_change, fill= anchor))+
#   geom_bar(stat = "count")+
#   ggtitle("difference between radius 100 - 50, if more cells at 100 -> postive value")+ # Why do we have more at 50
#   theme_classic()+
#   facet_wrap(.~ counted)

# KNN_ROI_overall_bivariate4 %>% 
#   filter(!is.na(degree_of_clustering_permutation)) %>% 
#   group_by(suid, anchor, counted) %>% 
#   mutate(n = n()) %>% 
#   filter(n == 2) %>% 
#   ggplot(aes(x= r, y= degree_of_clustering_permutation, color= suid, group = suid))+
#   geom_line()+
#   ggtitle("Ex: CD8 better radius at 50")+
#   theme_classic()+
#   facet_wrap(counted ~ anchor)+ 
#   theme(legend.position = "none")
# 
# KNN_ROI_overall_bivariate[23:24,] %>% 
#   ggplot(aes(x= r, y= degree_of_clustering_permutation), color= suid)+
#   geom_line()
# KNN_ROI_overall_bivariate[239:240,] %>% 
#   ggplot(aes(x= r, y= degree_of_clustering_permutation), color= suid)+
#   geom_line()
# 
# KNN_ROI_overall_bivariate %>% 
#   ggplot(aes(x= as.factor(r), y= degree_of_clustering_permutation))+
#   geom_boxplot()

########################################################################### II ### ICC
# Let's choose r = 50 for now
df_overall <- KNN_ROI_overall_bivariate3 %>% 
  filter(r == 50 & anchor != "DAPI (DAPI) Positive") %>% 
  select(-r) %>% 
  as.data.frame()
df_intra <- df_overall %>% 
  filter(annotation == "Intratumoral")
df_peri <- df_overall %>% 
  filter(annotation == "Peripheral")


library(psych)
vec_col <- c("theoretical_csr", "permuted_g")
ICCC_data <- data.frame(matrix(nrow = 1, ncol = 0))
ICC_data <- data.frame(matrix(nrow = 1, ncol = 0))
lb_data <- data.frame(matrix(nrow = 1, ncol = 0))
up_data <- data.frame(matrix(nrow = 1, ncol = 0))
fct_icc <- function(data) {
  for (i in vec_col) {
    raw_df <- data %>% select(suid, anchor, counted)
    
    if (class(data[, i]) == "numeric" |
        class(data[, i]) == "integer") {
      ICC_df <- cbind(raw_df, value = data[, i]) %>%
        drop_na() %>%
        filter(anchor == "CD3 (Opal 650) Positive" & counted == "CD3+ CD8+")
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
  colnames(ICCC_data) <- vec_col
  ICCC_data <- as.data.frame(t(ICCC_data)) %>%
    mutate(ICC_lb_up = 
             paste(round(V1, 2), " (", 
                   round(V2, 2), ", ", 
                   round(V3, 2), ")", sep = "")) %>%
    select(ICC_lb_up)
}
ICC_results <- fct_icc(df_overall)
ICC_results

ICC_results <- fct_icc(df_intra)
ICC_results

ICC_results <- fct_icc(df_peri)
ICC_results

rm(ICC_data, ICCC_data, lb_data, up_data, df, ICC_results)


########################################################################### III ### Summarize
# Because ICC is good, we can summarize slides measurement by suid, anchor, counted and annotation

KNN_ROI_overall_bivariate4 <- KNN_ROI_overall_bivariate3 %>% 
  filter(r == 50 & anchor != "DAPI (DAPI) Positive") %>% 
  group_by(suid, anchor, counted, annotation) %>% 
  summarize(theoretical_csr = mean(theoretical_csr, na.rm=TRUE),
            permuted_g = mean(permuted_g, na.rm=TRUE),
            observed_g = mean(observed_g, na.rm=TRUE),
            degree_of_clustering_permutation = mean(degree_of_clustering_permutation, na.rm=TRUE),
            degree_of_clustering_theoretical = mean(degree_of_clustering_theoretical, na.rm=TRUE)) %>% 
  ungroup() %>% 
  distinct(suid, anchor, counted, annotation, .keep_all = TRUE) %>% 
  mutate(across(where(is.numeric), ~ na_if(., NaN)))

KNN_ROI_overall_bivariate_intra <- KNN_ROI_overall_bivariate4 %>% 
  filter(annotation == "Intratumoral") %>% 
  select(-annotation)
KNN_ROI_overall_bivariate_peri <- KNN_ROI_overall_bivariate4 %>% 
  filter(annotation == "Peripheral") %>% 
  select(-annotation)

KNN_ROI_overall_bivariate4 <- 
  full_join(KNN_ROI_overall_bivariate_intra, 
            KNN_ROI_overall_bivariate_peri,
            by= c("suid", "anchor", "counted"),
            suffix = c("_i", "_p"))
########################################################################### IV ### Data exploratory
# Why do we see so little cells CD3, it should be more than CD11

is.na(KNN_ROI_overall_bivariate4$degree_of_clustering_permutation_i)

KNN_ROI_overall_bivariate4 %>% 
  ggplot(aes(x= counted, y= theoretical_csr_i, color= counted))+
  geom_boxplot()+
  theme_classic()+
  facet_wrap(. ~ anchor)+
  coord_flip()+ 
  theme(legend.position = "none")

KNN_ROI_overall_bivariate4 %>% 
  ggplot(aes(x= counted, y= degree_of_clustering_permutation_i, color= counted))+
  geom_violin()+
  theme_classic()+
  facet_wrap(. ~ anchor)+
  coord_flip()+ 
  theme(legend.position = "none")


########################################################################### V ### Make groups
summarized_markers_ROI <- readRDS("~/Documents/GitHub/Peres/spatial_contexture/summarized_markers_ROI.rds") %>% 
  mutate(vitalstatus = case_when(
    vitalstatus == "Deceased"                ~ 1,
    vitalstatus == "Alive"                   ~ 0
  )) %>% 
  select(suid, pair_id, 
         vitalstatus, timelastfu,
         refage, stage, 
         cd3_tumor_i : immunoscore_2018lancet_patients)

KNN_ROI_overall_bivariate5 <- KNN_ROI_overall_bivariate4 %>% 
  full_join(., summarized_markers_ROI, 
            by = "suid")

KNN_ROI_overall_bivariate6 <- KNN_ROI_overall_bivariate5 %>% 
  filter(anchor == "CD11b (Opal 620) Positive" & counted == "CD11b+ CD15+") %>% 
  mutate(CD11_vs_CD11CD15_grp = case_when(
    cd11b_total_i == "Absence"                            ~ "Absence",
    ntile(theoretical_csr_i, 2) == 1       ~ "Low",
    ntile(theoretical_csr_i, 2) == 2       ~ "High",
  ), CD11_vs_CD11CD15_grp = factor(CD11_vs_CD11CD15_grp, levels = c("Absence", "Low","High"))) %>% 
  select(suid, anchor, counted, CD11_vs_CD11CD15_grp)

KNN_ROI_overall_bivariate7 <- KNN_ROI_overall_bivariate5 %>% 
  filter(anchor == "CD3 (Opal 650) Positive" & counted == "CD3+ CD8+") %>% 
  mutate(tertile = ntile(theoretical_csr_i, 2)) %>% 
  mutate(CD3_vs_CD3CD8_grp = case_when(
    tertile == 1                                          ~ "Low",
    tertile == 2                                          ~ "High"
  ), CD3_vs_CD3CD8_grp = factor(CD3_vs_CD3CD8_grp, levels = c("Low","High"))) %>% 
  select(suid, anchor, counted, CD3_vs_CD3CD8_grp)

KNN_ROI_overall_bivariate_final <- KNN_ROI_overall_bivariate5 %>% 
  full_join(KNN_ROI_overall_bivariate6, ., 
            by= c("suid", "anchor", "counted")) %>% 
  full_join(KNN_ROI_overall_bivariate7, .,
            by= c("suid", "anchor", "counted"))


########################################################################### VI ### Survival

KNN_ROI_overall_bivariate_final

tbl1 <- KNN_ROI_overall_bivariate_final %>% 
  select(vitalstatus, timelastfu,
         CD11_vs_CD11CD15_grp,
         refage, race, stage) %>%
  tbl_uvregression(method = survival::coxph,
                   y = (Surv(time = KNN_ROI_overall_bivariate_final$timelastfu,
                             event = KNN_ROI_overall_bivariate_final$vitalstatus)),
                   exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>%
  bold_p(t = .05) %>% add_nevent(location = "level") %>% add_n(location = "level")
tbl2 <-
  coxph(Surv(time = KNN_ROI_overall_bivariate_final$timelastfu,
             event = KNN_ROI_overall_bivariate_final$vitalstatus) ~
          CD11_vs_CD11CD15_grp + refage + stage,
        data =  KNN_ROI_overall_bivariate_final) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_p(t = .05) %>%
  add_nevent(location = "level") %>% add_n(location = "level")
tbl3 <- KNN_ROI_overall_bivariate_final %>% filter(!is.na(CD11_vs_CD11CD15_grp))
tbl3 <- 
  coxph(Surv(time = tbl3$timelastfu,
             event = tbl3$vitalstatus) ~
          theoretical_csr_i + refage + stage,
        data =  tbl3) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_p(t = .05) %>%
  add_nevent(location = "level") %>% add_n(location = "level")
tbl_merge(list(tbl1, tbl2, tbl3), tab_spanner = c("**Univariable**", "**Multivariable**", "**Continuous**"))

tbl1 <- KNN_ROI_overall_bivariate_final %>% 
  select(vitalstatus, timelastfu,
         CD3_vs_CD3CD8_grp,
         refage, race, stage) %>%
  tbl_uvregression(method = survival::coxph,
                   y = (Surv(time = KNN_ROI_overall_bivariate_final$timelastfu,
                             event = KNN_ROI_overall_bivariate_final$vitalstatus)),
                   exponentiate = TRUE) %>%
  bold_labels() %>% italicize_levels() %>%
  bold_p(t = .05) %>% add_nevent(location = "level") %>% add_n(location = "level")
tbl2 <-
  coxph(Surv(time = KNN_ROI_overall_bivariate_final$timelastfu,
             event = KNN_ROI_overall_bivariate_final$vitalstatus) ~
          CD3_vs_CD3CD8_grp + refage + stage,
        data =  KNN_ROI_overall_bivariate_final) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_p(t = .05) %>%
  add_nevent(location = "level") %>% add_n(location = "level")
tbl3 <- KNN_ROI_overall_bivariate_final %>% filter(!is.na(CD3_vs_CD3CD8_grp))
tbl3 <- 
  coxph(Surv(time = tbl3$timelastfu,
             event = tbl3$vitalstatus) ~
          theoretical_csr_i + refage + stage,
        data =  tbl3) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_p(t = .05) %>%
  add_nevent(location = "level") %>% add_n(location = "level")
tbl_merge(list(tbl1, tbl2, tbl3), tab_spanner = c("**Univariable**", "**Multivariable**", "**Continuous**"))










