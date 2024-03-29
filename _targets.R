# Create path
path1 <- fs::path("","Volumes","Peres_Research", "K99_R00", "Image analysis data", "Spatial Data",
                 "Alex", "computed", "AACES_ROI1")
path2 <- fs::path("","Volumes","Peres_Research", "K99_R00", "Image analysis data", "Spatial Data",
                 "Alex", "computed", "AACES_ROI2")

# Import package
library(targets)

# Comfig
tar_option_set(packages = c("readr", "tidyverse")) # Package needed to run pipeline without importing them


# Run pipeline aka used to be the plan.R with drake
list(
  # ROI1
  tar_target(
    KNN_ROI1_overall_bivariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_overall_mif.rds")))[["derived"]][["bivariate_NN"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted G", "Observed G",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI1_overall_univariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_overall_mif.rds")))[["derived"]][["univariate_NN"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted CSR", "Observed",
             "Degree of Clustering Theoretical", "Degree of Clustering Permutation")
  ),
  tar_target(
    C_ROI1_overall_bivariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_overall_mif.rds")))[["derived"]][["bivariate_Count"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    C_ROI1_overall_univariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_overall_mif.rds")))[["derived"]][["univariate_Count"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI1_stroma_bivariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_stroma_mif.rds")))[["derived"]][["bivariate_NN"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted G", "Observed G",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI1_stroma_univariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_stroma_mif.rds")))[["derived"]][["univariate_NN"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted CSR", "Observed",
             "Degree of Clustering Theoretical", "Degree of Clustering Permutation")
  ),
  tar_target(
    C_ROI1_stroma_bivariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_stroma_mif.rds")))[["derived"]][["bivariate_Count"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    C_ROI1_stroma_univariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_stroma_mif.rds")))[["derived"]][["univariate_Count"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI1_tumor_bivariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_tumor_mif.rds")))[["derived"]][["bivariate_NN"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted G", "Observed G",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI1_tumor_univariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_tumor_mif.rds")))[["derived"]][["univariate_NN"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted CSR", "Observed",
             "Degree of Clustering Theoretical", "Degree of Clustering Permutation")
  ),
  tar_target(
    C_ROI1_tumor_bivariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_tumor_mif.rds")))[["derived"]][["bivariate_Count"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    C_ROI1_tumor_univariate,
    read_rds(paste0(path1,
                    c("/CD3_fixed_AACES_ROI1_tumor_mif.rds")))[["derived"]][["univariate_Count"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  # ROI2
  tar_target(
    KNN_ROI2_overall_bivariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_overall_mif.rds")))[["derived"]][["bivariate_NN"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted G", "Observed G",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI2_overall_univariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_overall_mif.rds")))[["derived"]][["univariate_NN"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted CSR", "Observed",
             "Degree of Clustering Theoretical", "Degree of Clustering Permutation")
  ),
  tar_target(
    C_ROI2_overall_bivariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_overall_mif.rds")))[["derived"]][["bivariate_Count"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    C_ROI2_overall_univariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_overall_mif.rds")))[["derived"]][["univariate_Count"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI2_stroma_bivariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_stroma_mif.rds")))[["derived"]][["bivariate_NN"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted G", "Observed G",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI2_stroma_univariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_stroma_mif.rds")))[["derived"]][["univariate_NN"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted CSR", "Observed",
             "Degree of Clustering Theoretical", "Degree of Clustering Permutation")
  ),
  tar_target(
    C_ROI2_stroma_bivariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_stroma_mif.rds")))[["derived"]][["bivariate_Count"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    C_ROI2_stroma_univariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_stroma_mif.rds")))[["derived"]][["univariate_Count"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI2_tumor_bivariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_tumor_mif.rds")))[["derived"]][["bivariate_NN"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted G", "Observed G",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    KNN_ROI2_tumor_univariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_tumor_mif.rds")))[["derived"]][["univariate_NN"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted CSR", "Observed",
             "Degree of Clustering Theoretical", "Degree of Clustering Permutation")
  ),
  tar_target(
    C_ROI2_tumor_bivariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_tumor_mif.rds")))[["derived"]][["bivariate_Count"]] %>%
      select("Image Location", "anchor", "counted", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  ),
  tar_target(
    C_ROI2_tumor_univariate,
    read_rds(paste0(path2,
                    c("/CD3_fixed_AACES_ROI2_tumor_mif.rds")))[["derived"]][["univariate_Count"]] %>%
      select("Image Location", "Marker", "r",
             "Theoretical CSR", "Permuted K", "Observed K",
             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  )
)


# tar_make() #### can NOT RUN in script, must be run in the console

# tar_load(KNN_ROI_overall_bivariate)
# tar_load(KNN_ROI_overall_univariate)
# tar_load(C_ROI_overall_bivariate)
# tar_load(C_ROI_overall_univariate)

# drake_clean =
# tar_destroy(destroy = "local")