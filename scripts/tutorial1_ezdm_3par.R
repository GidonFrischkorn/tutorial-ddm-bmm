#' Tutorial 1: EZ-Diffusion Model (3-Parameter) with bmm
#'
#' This tutorial demonstrates the complete workflow for fitting the 3-parameter
#' EZDM to aggregated summary statistics using bmm. The 3par version estimates
#' drift, bound, and ndt while fixing the relative starting point (zr) at 0.5.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

# -- Download pre-fitted model objects from OSF (optional) --------------------
download_from_osf <- FALSE
if (download_from_osf) source(here::here("scripts", "00_download_osf.R"))

pacman::p_load(
  here, bmm, brms, bridgesampling,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

# -- Load aggregated data ------------------------------------------------------
data_ezdm <- read_csv(here("data", "simulated_ezdm_3par.csv"))

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

# TODO: Specify the EZDM model
#
# model <- ezdm(
#   resp_vars = c("mean_rt", "var_rt", "n_upper", "n_trials"),
#   version = "3par"
# )
#
# # Define formula: condition effects on drift rate
# formula <- bmf(
#   drift ~ 0 + condition + (0 + condition | id),
#   bound ~ 1 + (1 | id),
#   ndt   ~ 1 + (1 | id)
# )
#
# # Set priors
# prior <- set_prior("normal(0, 1)", class = "b", dpar = "drift") +
#   set_prior("normal(0, 0.5)", class = "b", dpar = "bound") +
#   set_prior("normal(-1.5, 0.5)", class = "b", dpar = "ndt")

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

# TODO: Fit the model
#
# fit_ezdm_3par <- bmm(
#   formula = formula,
#   data = data_ezdm,
#   model = model,
#   prior = prior,
#   chains = chains,
#   cores = cores,
#   warmup = warmup,
#   iter = iter,
#   backend = "cmdstanr",
#   sample_prior = "yes",
#   file = here("output", "fit_ezdm_3par")
# )

###############################################################################!
# 3) Model Evaluation -----------------------------------------------------------
###############################################################################!

# TODO: Convergence diagnostics (Rhat, ESS, trace plots)
# TODO: Posterior predictive checks
# TODO: Parameter estimates and condition effects
# TODO: Hypothesis testing with brms::hypothesis()
