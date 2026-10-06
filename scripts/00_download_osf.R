#' Download Pre-Fitted Model Objects from OSF
#'
#' Run this script to download all fitted model objects (.rds) from OSF
#' into the output/ directory. This avoids re-fitting models from scratch,
#' which can take several hours.

pacman::p_load(here, osfr)

# Create output directory if it doesn't exist
dir.create(here("output"), showWarnings = FALSE)

# TODO: Set up OSF project and add download links
# osf_project <- osf_retrieve_node("XXXXX")
# osf_files <- osf_ls_files(osf_project, path = "output")
# osf_download(osf_files, path = here("output"), conflicts = "skip")

message("OSF download is not set up yet: nothing was downloaded. ",
        "Run the tutorial scripts to fit the models (see README).")
