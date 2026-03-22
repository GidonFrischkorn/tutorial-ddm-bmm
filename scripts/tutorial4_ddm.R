#' Tutorial 4: The Diffusion Decision Model and Cross-Model Comparison
#'
#' This tutorial fits the full DDM using the Wiener first-passage time
#' likelihood and then systematically compares parameter estimates across
#' all three EAM implementations (EZDM, CSWALD, DDM). Using simulated data
#' with known generating parameters, we evaluate how well each model
#' recovers the ground truth and compare computational cost. We also
#' demonstrate the DDM on a real published dataset.
#'
#' The central question: How do the three models compare, and when should
#' you use which?

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

# -- Load cleaned trial-level data (simulated) ---------------------------------
data_sim <- read_csv(here("data", "simulated_rt_data_clean.csv"))

# -- Load true generating parameters ------------------------------------------
true_params <- read_csv(here("data", "true_parameters.csv"))

# -- Load previous model fits for comparison -----------------------------------
# fit_ezdm_3par    <- readRDS(here("output", "fit_ezdm_3par.rds"))
# fit_cswald_crisk <- readRDS(here("output", "fit_cswald_crisk.rds"))

###############################################################################!
# 1) DDM Specification & Fitting -----------------------------------------------
###############################################################################!

# The DDM uses the Wiener first-passage time likelihood, which models the
# full RT distribution for both correct and error responses. This is the
# most complete EAM available in bmm but also the most computationally
# demanding.
#
# TODO: Specify and fit the DDM
#
# model <- ddm(
#   rt = "rt",
#   response = "response"
# )
#
# formula <- bmf(
#   drift ~ 0 + condition + (0 + condition | id),
#   bound ~ 1 + (1 | id),
#   ndt   ~ 1 + (1 | id)
# )
#
# t_start_ddm <- Sys.time()
# fit_ddm <- bmm(
#   formula = formula,
#   data = data_sim,
#   model = model,
#   prior = prior,
#   chains = chains, cores = cores,
#   warmup = warmup, iter = iter,
#   backend = "cmdstanr",
#   sample_prior = "yes",
#   file = here("output", "fit_ddm")
# )
# t_ddm <- difftime(Sys.time(), t_start_ddm, units = "mins")

## 1.1) Convergence diagnostics ------------------------------------------------
# TODO: Rhat, ESS, trace plots

## 1.2) Posterior predictive checks ---------------------------------------------
# TODO: pp_check for correct and error RT distributions
# The DDM should capture the full shape of both distributions, including
# the characteristic positive skew and the speed of errors vs. corrects.

###############################################################################!
# 2) Cross-Model Parameter Recovery --------------------------------------------
###############################################################################!

# This is the key comparison: all three models (EZDM, CSWALD, DDM) fitted to
# the same simulated data with known generating parameters. We compare how
# well each model recovers the true drift, bound, and ndt.
#
# TODO: Extract population-level estimates from all three models
#
# extract_estimates <- function(fit, model_name) {
#   fe <- fixef(fit)
#   # Extract and label relevant parameters
#   # (exact extraction depends on model parameterization)
#   tibble(
#     model = model_name,
#     parameter = rownames(fe),
#     estimate = fe[, "Estimate"],
#     lower = fe[, "Q2.5"],
#     upper = fe[, "Q97.5"]
#   )
# }
#
# all_estimates <- bind_rows(
#   extract_estimates(fit_ezdm_3par, "EZDM"),
#   extract_estimates(fit_cswald_crisk, "CSWALD"),
#   extract_estimates(fit_ddm, "DDM")
# )

## 2.1) Parameter recovery plot ------------------------------------------------
#
# Figure: Estimated vs. true parameter values for each model x parameter.
# Each panel shows one parameter (drift, bound, ndt). Within each panel,
# the three models are compared against the true generating value.
#
# TODO: Create parameter recovery figure
#
# p_recovery <- all_estimates |>
#   left_join(true_params_long, by = "parameter") |>
#   ggplot(aes(x = model, y = estimate, color = model)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
#   geom_hline(aes(yintercept = true_value),
#              linetype = "dashed", color = "grey50") +
#   facet_wrap(~parameter, scales = "free_y") +
#   scale_color_eam() +
#   labs(x = "Model", y = "Parameter Estimate",
#        caption = "Dashed line = true generating value") +
#   clean_plot()

## 2.2) Parameter recovery summary table ---------------------------------------
#
# Table: Bias (estimate - true), coverage (does 95% CI contain true value?),
# and precision (CI width) for each model x parameter combination.
#
# TODO: Compute recovery metrics
#
# recovery_table <- all_estimates |>
#   left_join(true_params_long, by = "parameter") |>
#   mutate(
#     bias = estimate - true_value,
#     coverage = lower <= true_value & true_value <= upper,
#     precision = upper - lower
#   )

###############################################################################!
# 3) Computational Cost Comparison ---------------------------------------------
###############################################################################!

# Compare fitting times across all three models. The CSWALD should be
# substantially faster than the DDM, while the EZDM (fitted to aggregated
# data) should be the fastest of all.
#
# TODO: Load fitting times from previous tutorials and compare
#
# # Load CSWALD fitting times (saved in Tutorial 3)
# # t_cswald_crisk <- ...
#
# timing_comparison <- tibble(
#   model = c("EZDM 3par", "CSWALD crisk", "DDM"),
#   fitting_time_min = c(t_ezdm, t_cswald_crisk, t_ddm),
#   data_type = c("Aggregated", "Trial-level", "Trial-level"),
#   n_parameters = c(3, 4, 3)
# )
#
# # Bar plot of fitting times
# p_timing <- timing_comparison |>
#   ggplot(aes(x = reorder(model, fitting_time_min),
#              y = as.numeric(fitting_time_min),
#              fill = model)) +
#   geom_col() +
#   scale_fill_eam() +
#   labs(x = "Model", y = "Fitting Time (minutes)") +
#   clean_plot()

###############################################################################!
# 4) Real Data Application -----------------------------------------------------
###############################################################################!

# Demonstrate the DDM on a real published dataset to show a practical
# application beyond simulated data.
#
# TODO: Identify and load real dataset (e.g., Wagenmakers et al., 2008)
#
# data_real <- read_csv(here("data", "..."))
#
# fit_ddm_real <- bmm(
#   formula = formula_real,
#   data = data_real,
#   model = ddm(rt = "rt", response = "response"),
#   ...
#   file = here("output", "fit_ddm_real")
# )
#
# # Full evaluation: convergence, pp_check, parameter interpretation

###############################################################################!
# 5) Practical Decision Tree ---------------------------------------------------
###############################################################################!

# Summary: When to use which model?
#
# | Criterion          | EZDM           | CSWALD          | DDM             |
# |--------------------|----------------|-----------------|-----------------|
# | Data type          | Aggregated     | Trial-level     | Trial-level     |
# | Speed              | Fastest        | Fast            | Slowest         |
# | Error rate needed  | Any            | <20% (simple),  | Any             |
# |                    |                | any (crisk)     |                 |
# | RT distribution    | No (summary    | Partial (correct| Full (correct + |
# | fit                | stats only)    | boundary only   | error)          |
# |                    |                | for simple)     |                 |
# | Best for           | Exploration,   | Large samples,  | Final analysis, |
# |                    | many conditions| iterative dev,  | publication,    |
# |                    |                | quick proxy     | full RT shape   |
# | Parameter recovery | Good (if data  | Good (proxy for | Gold standard   |
# |                    | quality is OK) | DDM parameters) |                 |
#
# Recommended workflow:
# 1. Start with EZDM for quick exploration and condition screening
# 2. Use CSWALD for iterative model development (fast turnaround)
# 3. Fit DDM for final analysis and publication
# 4. Compare CSWALD and DDM estimates as a robustness check
