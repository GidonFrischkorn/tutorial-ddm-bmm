#' Prepare Simulated Data for EAM Tutorials
#'
#' Generates trial-level data from known DDM parameters using rddm().
#' Produces two versions:
#'   1) Clean data — ground truth for parameter recovery
#'   2) Contaminated data — with injected fast guesses and slow outliers
#'
#' The contaminated data is used in Tutorial 0 to demonstrate how contaminants
#' bias EAM parameter estimates. All subsequent tutorials use the cleaned data.

###############################################################################!
# 0) R Setup: Packages ---------------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, tidyverse)

set.seed(2025)

###############################################################################!
# 1) Simulation Parameters -----------------------------------------------------
###############################################################################!

# Design: 50 participants, 2 conditions (easy vs hard), 200 trials per condition
n_subjects    <- 50
n_trials_cond <- 200
conditions    <- c("easy", "hard")

# Population-level DDM parameters
# Easy condition: higher drift rate -> faster, more accurate
# Hard condition: lower drift rate -> slower, less accurate
pop_params <- tibble(
  condition = conditions,
  drift     = c(2.5, 1.0),   # drift rate
  bound     = c(1.5, 1.5),   # boundary separation (same across conditions)
  ndt       = c(0.3, 0.3),   # non-decision time in seconds (same)
  zr        = c(0.5, 0.5)    # relative starting point (unbiased)
)

# Between-subject variability (SD of random effects)
tau_drift <- 0.5   # SD on identity scale
tau_bound <- 0.15  # SD on log scale
tau_ndt   <- 0.03  # SD on identity scale

###############################################################################!
# 2) Generate Clean Data -------------------------------------------------------
###############################################################################!

# TODO: Generate trial-level data using rddm() once bmm's DDM branch is merged.
# For now, this script provides the structure; actual data generation will use:
#
# data_clean <- expand_grid(
#   id = 1:n_subjects,
#   condition = conditions,
#   trial = 1:n_trials_cond
# ) |>
#   left_join(pop_params, by = "condition") |>
#   group_by(id) |>
#   mutate(
#     # Add subject-level random effects
#     drift_subj = drift + rnorm(1, 0, tau_drift),
#     bound_subj = bound * exp(rnorm(1, 0, tau_bound)),
#     ndt_subj   = ndt + rnorm(1, 0, tau_ndt)
#   ) |>
#   ungroup() |>
#   rowwise() |>
#   mutate(
#     sim = list(rddm(1, drift = drift_subj, bound = bound_subj,
#                      ndt = ndt_subj, zr = zr))
#   ) |>
#   unnest(sim) |>
#   select(id, condition, trial, rt, response,
#          drift_subj, bound_subj, ndt_subj)

###############################################################################!
# 3) Inject Contaminants -------------------------------------------------------
###############################################################################!

# Create a contaminated version of the data by replacing ~7% of trials with
# contaminant RTs. This mirrors realistic data where some responses arise from
# processes other than evidence accumulation (Ratcliff & Kang, 2021).
#
# Contaminant types:
#   a) Fast guesses (~4%): RT drawn from Uniform(0.1, ndt_subj), with random
#      accuracy at chance (Bernoulli 0.5). These represent responses made
#      before evidence accumulation could meaningfully begin.
#   b) Slow outliers (~3%): RT drawn from Uniform(mean_RT + 2*SD_RT,
#      mean_RT + 5*SD_RT). These represent attentional lapses or delayed
#      processing unrelated to the decision.
#
# TODO: Implement contaminant injection once clean data generation works.
#
# contam_rate_fast <- 0.04
# contam_rate_slow <- 0.03
#
# data_contaminated <- data_clean |>
#   group_by(id, condition) |>
#   mutate(
#     n = n(),
#     is_fast_contam = row_number() <= round(n * contam_rate_fast),
#     is_slow_contam = row_number() > (n - round(n * contam_rate_slow)),
#     # Shuffle which trials are contaminated
#     contam_order = sample(n)
#   ) |>
#   arrange(id, condition, contam_order) |>
#   mutate(
#     is_fast_contam = row_number() <= round(n * contam_rate_fast),
#     is_slow_contam = !is_fast_contam &
#       row_number() > (n - round(n * contam_rate_slow))
#   ) |>
#   ungroup() |>
#   mutate(
#     mean_rt_group = mean(rt),
#     sd_rt_group   = sd(rt),
#     # Inject fast guesses
#     rt = if_else(is_fast_contam,
#                  runif(n(), 0.1, ndt_subj),
#                  rt),
#     response = if_else(is_fast_contam,
#                        as.numeric(rbinom(n(), 1, 0.5)),
#                        response),
#     # Inject slow outliers
#     rt = if_else(is_slow_contam,
#                  runif(n(), mean_rt_group + 2 * sd_rt_group,
#                        mean_rt_group + 5 * sd_rt_group),
#                  rt),
#     is_contaminant = is_fast_contam | is_slow_contam
#   ) |>
#   select(id, condition, trial, rt, response, is_contaminant)

###############################################################################!
# 4) Save Data -----------------------------------------------------------------
###############################################################################!

# TODO: Save all data versions once generation is implemented.
#
# # Clean trial-level data (ground truth)
# write_csv(
#   data_clean |> select(id, condition, trial, rt, response),
#   here("data", "simulated_rt_data_clean.csv")
# )
#
# # Contaminated trial-level data (for Tutorial 0 bias demonstration)
# write_csv(
#   data_contaminated,
#   here("data", "simulated_rt_data_contaminated.csv")
# )
#
# # Aggregated summary statistics for EZDM (3par) — from clean data
# data_ezdm_3par <- ezdm_summary_stats(
#   data_clean, rt = rt, response = response,
#   .by = c(id, condition), version = "3par"
# )
# write_csv(data_ezdm_3par, here("data", "simulated_ezdm_3par.csv"))
#
# # Aggregated summary statistics for EZDM (4par) — from clean data
# data_ezdm_4par <- ezdm_summary_stats(
#   data_clean, rt = rt, response = response,
#   .by = c(id, condition), version = "4par"
# )
# write_csv(data_ezdm_4par, here("data", "simulated_ezdm_4par.csv"))
#
# # Save true generating parameters for parameter recovery comparisons
# write_csv(pop_params, here("data", "true_parameters.csv"))
