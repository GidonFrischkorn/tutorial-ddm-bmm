# Download the fitted model objects from OSF into output/, so that the tutorial
# scripts load them through bmm()'s file argument instead of re-fitting
# (the trial-level fits take 3 to 12 hours each).

pacman::p_load(here, osfr)

# A personal access token is only needed for private components; set it in the
# environment (OSF_PAT), never in this file
if (nzchar(Sys.getenv("OSF_PAT"))) osf_auth(Sys.getenv("OSF_PAT"))

osf_project <- osf_retrieve_node("XXXXX")

dir.create(here("output"), showWarnings = FALSE)
dir.create(here("figures"), showWarnings = FALSE)

osf_download(osf_ls_files(osf_project, path = "output", n_max = Inf),
             path = here("output"), conflicts = "skip", progress = TRUE)
osf_download(osf_ls_files(osf_project, path = "figures", n_max = Inf),
             path = here("figures"), conflicts = "skip", progress = TRUE)
