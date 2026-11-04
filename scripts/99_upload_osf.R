# Upload the fitted model objects, figures, and rendered manuscript to the OSF
# project (https://osf.io/7y54v). Maintainer script: readers use
# 00_download_osf.R instead.
#
# Uploads only the current files: the top-level .rds in output/, the top-level
# .pdf in figures/, and the rendered manuscript and supplement PDFs. Archive
# subfolders (output/_archive_*, figures/_archive_*) stay local.
#
# Needs an OSF personal access token with write scope in the environment
# (OSF_PAT, e.g. in ~/.Renviron), never in this file.

pacman::p_load(here, osfr)

pat <- Sys.getenv("OSF_PAT")
if (!nzchar(pat)) stop("OSF_PAT is not set; add a write-scope token to ~/.Renviron")
osf_auth(pat)

osf_project <- osf_retrieve_node("7y54v")

uploads <- list(
  output     = list.files(here("output"), pattern = "\\.rds$", full.names = TRUE),
  figures    = list.files(here("figures"), pattern = "\\.pdf$", full.names = TRUE),
  manuscript = here("manuscript", c("tutorial-ddm-bmm.pdf", "supplement.pdf"))
)

stopifnot(all(file.exists(unlist(uploads))))

for (folder in names(uploads)) {
  files <- uploads[[folder]]
  message(sprintf("Uploading %d files to %s/", length(files), folder))
  osf_dir <- osf_mkdir(osf_project, path = folder)
  osf_upload(osf_dir, path = files, conflicts = "overwrite", progress = TRUE)
}

# Check: every local file is on OSF with the same size
remote <- purrr::map_dfr(names(uploads), function(folder) {
  ls <- osf_ls_files(osf_project, path = folder, n_max = Inf)
  data.frame(folder = folder, name = ls$name,
             remote_size = purrr::map_dbl(ls$meta, \(m) m$attributes$size))
})
local <- data.frame(folder = rep(names(uploads), lengths(uploads)),
                    name = basename(unlist(uploads)),
                    local_size = file.size(unlist(uploads)))
check <- merge(local, remote, all.x = TRUE)
print(check, row.names = FALSE)
n_ok <- sum(check$local_size == check$remote_size, na.rm = TRUE)
message(sprintf("UPLOAD CHECK: %d of %d files match in size", n_ok, nrow(check)))
