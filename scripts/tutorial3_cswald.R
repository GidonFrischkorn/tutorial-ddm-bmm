#' Tutorial 3: The Censored-Shifted Wald Model as a Fast DDM Proxy
#'
#' The CSWALD provides a computationally efficient alternative to the full
#' Diffusion Decision Model (DDM). By modeling the correct-response boundary
#' with a shifted Wald distribution and treating errors as either censored
#' observations (simple version) or a competing accumulator (crisk version),
#' the CSWALD fits substantially faster than the DDM while recovering similar
#' parameter estimates when used appropriately.
#'
#' This tutorial demonstrates both CSWALD versions on the same trial-level
#' data used across all tutorials, and compares estimates against the EZDM
#' results from Tutorial 1.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

# -- Load cleaned trial-level data ---------------------------------------------
data_rt <- read_csv(here("data", "simulated_rt_data_clean.csv"))

# -- Load true generating parameters for comparison ----------------------------
true_params <- read_csv(here("data", "true_parameters.csv"))

###############################################################################!
# 1) Model Specification: Simple Version ---------------------------------------
###############################################################################!

# The simple CSWALD treats errors as censored correct responses — i.e., trials
# where the accumulator had not yet reached the correct boundary when the
# (unobserved) error boundary was hit. This is appropriate when error rates
# are low (<20%), because the censoring assumption becomes less tenable as
# errors increase.
#
# Note: In the simple version, `bound` represents the distance from the
# starting point to the correct boundary (half of the full boundary
# separation in DDM terms).
#
# TODO: Specify the simple CSWALD model
#
# model_simple <- cswald(
#   rt = "rt",
#   response = "response",
#   version = "simple"
# )
#
# formula_simple <- bmf(
#   drift ~ 0 + condition + (0 + condition | id),
#   bound ~ 1 + (1 | id),
#   ndt   ~ 1 + (1 | id)
# )
#
# prior_simple <- set_prior("normal(0, 1)", class = "b", dpar = "drift") +
#   set_prior("normal(0, 0.5)", class = "b", dpar = "bound") +
#   set_prior("normal(-2, 0.3)", class = "b", dpar = "ndt")

###############################################################################!
# 2) Model Fitting: Simple Version ---------------------------------------------
###############################################################################!

# TODO: Fit the simple CSWALD model and record fitting time
#
# t_start_simple <- Sys.time()
# fit_cswald_simple <- bmm(
#   formula = formula_simple,
#   data = data_rt,
#   model = model_simple,
#   prior = prior_simple,
#   chains = chains, cores = cores,
#   warmup = warmup, iter = iter,
#   backend = "cmdstanr",
#   sample_prior = "yes",
#   file = here("output", "fit_cswald_simple")
# )
# t_cswald_simple <- difftime(Sys.time(), t_start_simple, units = "mins")
# message("CSWALD simple fitting time: ", round(t_cswald_simple, 1), " minutes")

###############################################################################!
# 3) Model Specification & Fitting: Competing-Risks Version --------------------
###############################################################################!

# The competing-risks (crisk) version models both correct and error responses
# as separate racing accumulators. This is more appropriate when error rates
# are substantial (>10-20%), as it explicitly models the error process rather
# than treating it as censored.
#
# In the crisk version, `bound` represents the full boundary separation
# (consistent with DDM), and `zr` captures starting-point bias.
#
# TODO: Specify and fit the crisk CSWALD model
#
# model_crisk <- cswald(
#   rt = "rt",
#   response = "response",
#   version = "crisk"
# )
#
# formula_crisk <- bmf(
#   drift ~ 0 + condition + (0 + condition | id),
#   bound ~ 1 + (1 | id),
#   ndt   ~ 1 + (1 | id),
#   zr    ~ 1 + (1 | id)
# )
#
# t_start_crisk <- Sys.time()
# fit_cswald_crisk <- bmm(
#   formula = formula_crisk,
#   data = data_rt,
#   model = model_crisk,
#   prior = prior_crisk,
#   chains = chains, cores = cores,
#   warmup = warmup, iter = iter,
#   backend = "cmdstanr",
#   sample_prior = "yes",
#   file = here("output", "fit_cswald_crisk")
# )
# t_cswald_crisk <- difftime(Sys.time(), t_start_crisk, units = "mins")
# message("CSWALD crisk fitting time: ", round(t_cswald_crisk, 1), " minutes")

###############################################################################!
# 4) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 4.1) Convergence diagnostics ------------------------------------------------
#
# TODO: Check Rhat, ESS, trace plots for both models

## 4.2) Posterior predictive checks on RT distributions -------------------------
#
# TODO: pp_check for both models — compare predicted vs observed RT densities
#       for correct and error responses separately

## 4.3) Compare simple vs crisk parameter estimates ----------------------------
#
# TODO: Extract and compare population-level estimates
# - drift estimates should be similar across versions
# - bound differs in interpretation (half vs full boundary)
# - crisk additionally estimates zr

## 4.4) Compare CSWALD vs EZDM estimates ---------------------------------------
#
# Load the EZDM fit from Tutorial 1 and compare parameter estimates on the
# same data. This sets up the three-way comparison completed in Tutorial 4.
#
# TODO: Load EZDM fit and compare
#
# fit_ezdm_3par <- readRDS(here("output", "fit_ezdm_3par.rds"))
#
# # Extract population-level estimates from both models
# ezdm_estimates <- fixef(fit_ezdm_3par)
# cswald_estimates <- fixef(fit_cswald_crisk)
#
# # Compare: Do the CSWALD and EZDM agree on drift, bound, ndt?
# # They should be similar, demonstrating that the simpler EZDM captures the
# # same effects as the trial-level CSWALD — but the CSWALD additionally
# # provides full posterior predictive checks on the RT distribution shape.

## 4.5) Computational cost summary ---------------------------------------------
#
# Record fitting times. The CSWALD should fit substantially faster than the
# DDM (Tutorial 4) while providing similar estimates. This is the key
# advantage: when full DDM fitting is too slow (e.g., many conditions,
# large samples, iterative model development), the CSWALD serves as an
# efficient proxy.
#
# TODO: Print fitting time comparison
# message("Fitting times:")
# message("  CSWALD simple: ", round(t_cswald_simple, 1), " minutes")
# message("  CSWALD crisk:  ", round(t_cswald_crisk, 1), " minutes")
# message("  (DDM times will be compared in Tutorial 4)")
