# Download pre-fitted model objects from OSF
# Saves all fitted models (.rds) to output/, so the tutorial scripts load them
# instead of re-fitting (which takes several hours).

pacman::p_load(here, osfr)

dir.create(here("output"), showWarnings = FALSE)

# osf_project <- osf_retrieve_node("XXXXX")
# osf_files <- osf_ls_files(osf_project, path = "output")
# osf_download(osf_files, path = here("output"), conflicts = "skip")

message("OSF download is not set up yet: nothing was downloaded. ",
        "Run the tutorial scripts to fit the models (see README).")
