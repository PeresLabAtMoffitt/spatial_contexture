# Create path
path <- fs::path("","Volumes","Peres_Research", "K99_R00", "Image analysis data", "Spatial Data")

# Load package
library(targets)

# Comfig
tar_option_set(packages = c("readr", "tidyverse")) # Package needed to run pipeline without importing them


# Run pipeline aka used to be the plan.R with drake
list(
  # tar_target(
  #   raw_data_file,
  #   paste0(path,
  #   "/op.csv"),
  #   format = "file"
  # ),
  tar_target(
    ROI_overall_bivariate,
    read_csv(paste0(path,
                    "/Nearest Neighbor/ROI_overall_bivariate.csv"), col_types = cols()) %>%#________________
      select(suid, slide_type, image.tag, Marker : ncol(.))
  ),
  tar_target(
    ROI_stroma_bivariate,
    read_csv(paste0(path,
                    "/Nearest Neighbor/ROI_stroma_bivariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, anchor : ncol(.))
  ),
  tar_target(
    ROI_stroma_univariate,
    read_csv(paste0(path,
                    "/Nearest Neighbor/ROI_stroma_univariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, Marker : ncol(.))
  ),
  tar_target(
    ROI_tumor_bivariate,
    read_csv(paste0(path,
                    "/Nearest Neighbor/ROI_tumor_bivariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, anchor : ncol(.))
  ),
  tar_target(
    ROI_tumor_univariate,
    read_csv(paste0(path,
                    "/Nearest Neighbor/ROI_tumor_univariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, Marker : ncol(.))
  ),
  tar_target(
    K_ROI_overall_bivariate,
    read_csv(paste0(path,
                    "/Ripley_K/ROI_overall_bivariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, Marker : ncol(.))
  ),
  tar_target(
    K_ROI_overall_univariate,
    read_csv(paste0(path,
                    "/Ripley_K/ROI_overall_univariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, Marker : ncol(.))
  ),
  tar_target(
    K_ROI_stroma_bivariate,
    read_csv(paste0(path,
                    "/Ripley_K/ROI_stroma_bivariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, anchor : ncol(.))
  ),
  tar_target(
    K_ROI_stroma_univariate,
    read_csv(paste0(path,
                    "/Ripley_K/ROI_stroma_univariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, Marker : ncol(.))
  ),
  tar_target(
    K_ROI_tumor_bivariate,
    read_csv(paste0(path,
                    "/Ripley_K/ROI_tumor_bivariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, anchor : ncol(.))
  ),
  tar_target(
    K_ROI_tumor_univariate,
    read_csv(paste0(path,
                    "/Ripley_K/ROI_tumor_univariate.csv"), col_types = cols()) %>%
      select(suid, slide_type, image.tag, Marker : ncol(.))
  )
)


# tar_make() #### can NOT RUN in script, must be run in the console

# tar_load(ROI_overall_bivariate) #### can NOT RUN in script, must be run in the console
# tar_load(ROI_stroma_bivariate)
# tar_load(ROI_stroma_univariate)
# tar_load(ROI_tumor_bivariate)
# tar_load(ROI_tumor_univariate)
# 
# tar_load(K_ROI_overall_bivariate)
# tar_load(K_ROI_overall_univariate)
# tar_load(K_ROI_stroma_bivariate)
# tar_load(K_ROI_stroma_univariate)
# tar_load(K_ROI_tumor_bivariate)
# tar_load(K_ROI_tumor_univariate)

# drake_clean =
# tar_destroy()