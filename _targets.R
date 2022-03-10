# Create path
path <- fs::path("","Volumes","Peres_Research", "K99_R00", "Image analysis data")

# Load package
library(targets)

# Comfig
tar_option_set(packages = c("readr")) # Package needed to run pipeline without importing them


# Run pipeline aka used to be the plan.R with drake
list(
  tar_target(
    raw_data_file,
    paste0(path,
    "/op.csv"),
    format = "file"
  ),
  tar_target(
    raw_data,
    read_csv(raw_data_file, col_types = cols())
  )
)


# tar_make() #### NOT RUN in script but must be run in the console
# tar_load(raw_data) #### NOT RUN in script but must be run in the console

# drake_clean =
# tar_destroy()