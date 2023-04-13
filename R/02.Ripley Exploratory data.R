# Import package
library(targets)
library(tidyverse)
library(gtsummary)
library(survival)
library(survminer)

# Load data
tar_load(KNN_ROI_overall_univariate)
# tar_load(C_ROI_overall_bivariate)
# tar_load(C_ROI_overall_univariate)

KNN_ROI_overall_bivariate1 <- KNN_ROI_overall_univariate %>% 
  janitor::clean_names() %>% 
  ungroup() %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))

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


########################################################################### I ### AUC based on degree of clustering
# area under the Degree of Clustering curve (computed at various r) for choosing the best radius (r)
table(KNN_ROI_overall_bivariate3$r)


########################################################################### II ### ICC
# Let's choose r = 50 for now
df_overall <- KNN_ROI_overall_bivariate3 %>% 
  filter(r == 50 & marker != "DAPI (DAPI) Positive") %>% 
  select(-r) %>% 
  as.data.frame()
df_intra <- df_overall %>% 
  filter(annotation == "Intratumoral")
df_peri <- df_overall %>% 
  filter(annotation == "Peripheral")


library(psych)
vec_col <- c("theoretical_csr", "permuted_csr", "observed",
             "degree_of_clustering_theoretical", "degree_of_clustering_permutation")
# ICCC_data <- data.frame(matrix(nrow = 1, ncol = 0))
ICC_data <- data.frame(matrix(nrow = 1, ncol = 0))
lb_data <- data.frame(matrix(nrow = 1, ncol = 0))
up_data <- data.frame(matrix(nrow = 1, ncol = 0))
fct_icc <- function(data, cell_type) {
  for (i in vec_col) {
    raw_df <- data %>% select(suid, marker)
    
    if (class(data[, i]) == "numeric" |
        class(data[, i]) == "integer") {
      ICC_df <- cbind(raw_df, value = data[, i]) %>%
        drop_na() %>%
        filter(marker == cell_type)
      ICC_df <- ICC_df %>%
        mutate(Slide = "Slide0") %>%
        group_by(suid) %>%
        mutate(n = row_number(suid)) %>%
        ungroup() %>%
        unite(slide_id, Slide:n, sep = "", remove = TRUE, na.rm = TRUE) %>%
        pivot_wider(id_cols = -c(marker),
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
# ICC_results <- fct_icc(df_overall, cell_type = "CD11b (Opal 620) Positive")
# ICC_results
# 
# ICC_results <- fct_icc(df_intra, cell_type = "CD11b (Opal 620) Positive")
# ICC_results
# 
# ICC_results <- fct_icc(df_peri, cell_type = "CD11b (Opal 620) Positive")
# ICC_results

vec_mar <- c(unique(df_overall$marker))
ICC_lb_up_df_a <- data.frame(matrix(nrow = 1, ncol = 0))
ICC_lb_up_df_b <- data.frame(matrix(nrow = 1, ncol = 0))
ICC_lb_up_df_c <- data.frame(matrix(nrow = 1, ncol = 0))
for (i in vec_mar) {
  a <- fct_icc(df_overall, cell_type = i) %>% 
    `colnames<-`(i)
  ICC_lb_up_df_a <- cbind(ICC_lb_up_df_a, a) %>% 
    mutate(annotation = "Overall", .before= 1)
  b <- fct_icc(df_intra, cell_type = i) %>% 
    `colnames<-`(i)
  ICC_lb_up_df_b <- cbind(ICC_lb_up_df_b, b) %>% 
    mutate(annotation = "Intratumoral")
  c <- fct_icc(df_peri, cell_type = i) %>% 
    `colnames<-`(i)
  ICC_lb_up_df_c <- cbind(ICC_lb_up_df_c, c) %>% 
    mutate(annotation = "Peripheral")
  
}
bind_rows(ICC_lb_up_df_a,
          ICC_lb_up_df_b,
          ICC_lb_up_df_c)

rm(ICC_data, lb_data, up_data, a, b, c,
   ICC_lb_up_df_a, ICC_lb_up_df_b, ICC_lb_up_df_c,
   df_overall, df_intra, df_peri)


########################################################################### III ### Summarize
# Because ICC is good, we can summarize slides measurement by suid, anchor, counted and annotation

KNN_ROI_overall_bivariate4 <- KNN_ROI_overall_bivariate3 %>% 
  filter(r == 50) %>% 
  group_by(suid, marker, annotation) %>% 
  summarize(theoretical_csr = mean(theoretical_csr, na.rm=TRUE),
            permuted_csr = mean(permuted_csr, na.rm=TRUE),
            observed = mean(observed, na.rm=TRUE),
            degree_of_clustering_permutation = mean(degree_of_clustering_permutation, na.rm=TRUE),
            degree_of_clustering_theoretical = mean(degree_of_clustering_theoretical, na.rm=TRUE)) %>% 
  ungroup() %>% 
  distinct(suid, marker, annotation, .keep_all = TRUE) %>% 
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
            by= c("suid", "marker"),
            suffix = c("_i", "_p"))


########################################################################### IV ### Data exploratory
# Why do we see so little cells CD3, it should be more than CD11

is.na(KNN_ROI_overall_bivariate4$degree_of_clustering_permutation_i)

KNN_ROI_overall_bivariate4 %>% 
  ggplot(aes(x= marker, y= theoretical_csr_i, color= marker))+
  geom_boxplot()+
  theme_classic()+
  coord_flip()+ 
  theme(legend.position = "none")

KNN_ROI_overall_bivariate4 %>% 
  ggplot(aes(x= marker, y= degree_of_clustering_permutation_i, color= marker))+
  geom_violin()+
  theme_classic()+
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
         refage, race, stage, 
         cd3_tumor_i : immunoscore_2018lancet_patients)

KNN_ROI_overall_bivariate5 <- KNN_ROI_overall_bivariate4 %>% 
  full_join(., summarized_markers_ROI, 
            by = "suid")

KNN_ROI_overall_bivariate6 <- KNN_ROI_overall_bivariate5 %>% 
  filter(marker == "CD11b (Opal 620) Positive") %>% 
  mutate(spat_CD11_i = case_when(
    cd11b_total_i == "Absence"                            ~ "Absence",
    ntile(theoretical_csr_i, 2) == 1       ~ "Low",
    ntile(theoretical_csr_i, 2) == 2       ~ "High",
  ), spat_CD11_i = factor(spat_CD11_i, levels = c("Absence", "Low","High"))) %>% 
  select(suid, marker, spat_CD11_i)

KNN_ROI_overall_bivariate7 <- KNN_ROI_overall_bivariate5 %>% 
  filter(marker == "CD3 (Opal 650) Positive") %>% 
  mutate(tertile = ntile(theoretical_csr_i, 2)) %>% 
  mutate(spat_CD3_i = case_when(
    tertile == 1                                          ~ "Low",
    tertile == 2                                          ~ "High"
  ), spat_CD3_i = factor(spat_CD3_i, levels = c("Low","High"))) %>% 
  select(suid, marker, spat_CD3_i)

KNN_ROI_overall_bivariate_final <- KNN_ROI_overall_bivariate5 %>% 
  full_join(KNN_ROI_overall_bivariate6, ., 
            by= c("suid", "marker")) %>% 
  full_join(KNN_ROI_overall_bivariate7, .,
            by= c("suid", "marker"))

library(ggridges)
KNN_ROI_overall_bivariate_final %>% 
  select(spat_CD3_i, marker, theoretical_csr_i) %>% 
  # filter(!is.na(clusters_all_IandP)) %>% 
  # select(clusters_all_IandP, 
  #        percent_CD3 = percent_CD3_tumor.i, 
  #        percent_CD3_CD8 = percent_CD3_CD8_tumor.i, 
  #        percent_CD3_FoxP3 = percent_CD3_FoxP3_tumor.i, 
  #        percent_CD11b = percent_CD11b_tumor.i, 
  #        percent_CD11b_CD15 = percent_CD11b_CD15_tumor.i) %>% 
  # pivot_longer(-clusters_all_IandP, names_to = "markers", values_to = "value") %>% 
  mutate(location = "Tumor") %>% 
  # bind_rows(., b) %>% 
  
  ggplot(aes(y=spat_CD3_i, x=theoretical_csr_i,  fill=marker#, linetype = location
             )) +
  geom_density_ridges(alpha=0.5) +
  ggtitle("Intratumoral")+
  theme_ridges()+ facet_wrap(.~ marker)
  # theme(axis.text.y = element_blank(),
  #       panel.spacing = unit(0.1, "lines"),
  #       strip.text.x = element_text(size = 8)
  # )


########################################################################### VI ### Survival
tbl1 <- KNN_ROI_overall_bivariate_final %>% 
  select(vitalstatus, timelastfu,
         spat_CD11_i,
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
          spat_CD11_i + refage + stage,
        data =  KNN_ROI_overall_bivariate_final) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_p(t = .05) %>%
  add_nevent(location = "level") %>% add_n(location = "level")
# tbl3 <- KNN_ROI_overall_bivariate_final %>% filter(!is.na(spat_CD11_i))
# tbl3 <- 
#   coxph(Surv(time = tbl3$timelastfu,
#              event = tbl3$vitalstatus) ~
#           theoretical_csr_i + refage + stage,
#         data =  tbl3) %>%
#   tbl_regression(exponentiate = TRUE) %>%
#   bold_p(t = .05) %>%
#   add_nevent(location = "level") %>% add_n(location = "level")
tbl_merge(list(tbl1, tbl2), tab_spanner = c("**Univariable**", "**Multivariable**"))

tbl1 <- KNN_ROI_overall_bivariate_final %>% 
  select(vitalstatus, timelastfu,
         spat_CD3_i,
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
          spat_CD3_i + refage + stage,
        data =  KNN_ROI_overall_bivariate_final) %>%
  tbl_regression(exponentiate = TRUE) %>%
  bold_p(t = .05) %>%
  add_nevent(location = "level") %>% add_n(location = "level")
# tbl3 <- KNN_ROI_overall_bivariate_final %>% filter(!is.na(spat_CD3_i))
# tbl3 <- 
#   coxph(Surv(time = tbl3$timelastfu,
#              event = tbl3$vitalstatus) ~
#           theoretical_csr_i + refage + stage,
#         data =  tbl3) %>%
#   tbl_regression(exponentiate = TRUE) %>%
#   bold_p(t = .05) %>%
#   add_nevent(location = "level") %>% add_n(location = "level")
tbl_merge(list(tbl1, tbl2), tab_spanner = c("**Univariable**", "**Multivariable**"))

















