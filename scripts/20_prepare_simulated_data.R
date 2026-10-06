#' Prepare Simulated Data for EAM Tutorials
#'
#' Generates trial-level data from known DDM parameters using rddm().
#' Produces two datasets:
#'   1) An easy/hard design with an unbiased starting point -- ground truth
#'      for the parameter recovery study in the appendix (script 21)
#'   2) A stimulus x bias design with a biased starting point (zr = 0.60 or
#'      0.40) -- ground truth for the recovery study in Part I (script 22)
#'
#' The second dataset is generated after the first, so adding it left the
#' first dataset's draws unchanged (verified by diffing the written CSVs).

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

# Generate subject-level parameters by adding random effects to population values
subject_params <- expand_grid(
  id = 1:n_subjects,
  condition = conditions
) |>
  left_join(pop_params, by = "condition") |>
  group_by(id) |>
  mutate(
    drift_subj = drift + rnorm(1, 0, tau_drift),
    bound_subj = bound * exp(rnorm(1, 0, tau_bound)),
    ndt_subj   = pmax(0.1, ndt + rnorm(1, 0, tau_ndt))
  ) |>
  ungroup()

# Generate trial-level data using rddm()
# rowwise() is needed because each row (participant x condition) carries its
# own parameter values; reframe() then expands each row into n_trials_cond
# simulated trials
data_clean <- subject_params |>
  rowwise() |>
  reframe(
    id = id,
    condition = condition,
    trial = 1:n_trials_cond,
    rddm(n_trials_cond,
          drift = drift_subj, bound = bound_subj,
          ndt = ndt_subj, zr = zr),
    drift_true = drift_subj,
    bound_true = bound_subj,
    ndt_true   = ndt_subj
  )


###############################################################################!
# 3) Compute EZDM Summary Statistics -------------------------------------------
###############################################################################!

# Compute summary statistics from the clean data for EZDM tutorials

# 3-parameter version (drift, bound, ndt; zr fixed at 0.5)
data_ezdm_3par <- data_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par", method = "simple"),
    .by = c(id, condition)
  )

# 4-parameter version (adds zr via separate upper/lower moments)
data_ezdm_4par <- data_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par", method = "simple"),
    .by = c(id, condition)
  )

###############################################################################!
# 4) Save Data -----------------------------------------------------------------
###############################################################################!

# Clean trial-level data (ground truth for parameter recovery)
write_csv(
  data_clean |> select(id, condition, trial, rt, response),
  here("data", "simulated_rt_data_clean.csv")
)

# Aggregated summary statistics for EZDM (3par) -- from clean data
write_csv(data_ezdm_3par, here("data", "simulated_ezdm_3par.csv"))

# Aggregated summary statistics for EZDM (4par) -- from clean data
write_csv(data_ezdm_4par, here("data", "simulated_ezdm_4par.csv"))

# Save true generating parameters for parameter recovery comparisons
write_csv(pop_params, here("data", "true_parameters.csv"))

# Save subject-level true parameters for individual-level recovery
write_csv(
  subject_params |>
    select(id, condition, drift_subj, bound_subj, ndt_subj),
  here("data", "true_parameters_subject.csv")
)

###############################################################################!
# 5) Simulated Data With a Biased Starting Point -------------------------------
###############################################################################!

# Second dataset, for the recovery study in Part I of the manuscript. Two
# stimulus types set the direction of drift (+delta for A, -delta for B), so
# that with response coding (upper = an "A" response) both boundaries are
# reached often; two bias conditions set the relative starting point toward
# the "A" boundary (zr = 0.60) or the "B" boundary (zr = 0.40). Bound and ndt
# are those of the design above and do not vary by condition.
n_trials_bias <- 250

bias_params <- expand_grid(
  stimulus = c("A", "B"),
  bias     = c("A", "B")
) |>
  mutate(
    drift = if_else(stimulus == "A", 0.6, -0.6),
    bound = 1.5,
    ndt   = 0.3,
    zr    = if_else(bias == "A", 0.60, 0.40),
    n_trials_cell = n_trials_bias
  )

# Between-subject spreads: 0.25 on |drift|, 0.15 on log bound, 0.03 s on ndt,
# and 0.2 on the logit of zr. Drift varies less than in the design above so
# that every cell keeps enough responses at its minority boundary for the
# boundary-specific moments the 4par EZ-Diffusion Model needs.
tau_drift_bias <- 0.25
tau_zr_bias    <- 0.20

subject_bias <- tibble(id = 1:n_subjects) |>
  mutate(
    drift_re   = rnorm(n(), 0, tau_drift_bias),
    bound_subj = 1.5 * exp(rnorm(n(), 0, tau_bound)),
    ndt_subj   = pmax(0.1, 0.3 + rnorm(n(), 0, tau_ndt)),
    zr_re      = rnorm(n(), 0, tau_zr_bias)
  ) |>
  expand_grid(bias_params) |>
  mutate(
    drift_subj = sign(drift) * (abs(drift) + drift_re),
    zr_subj    = plogis(qlogis(zr) + zr_re)
  )

data_bias <- subject_bias |>
  rowwise() |>
  reframe(
    id = id, stimulus = stimulus, bias = bias,
    trial = 1:n_trials_bias,
    rddm(n_trials_bias, drift = drift_subj, bound = bound_subj,
         ndt = ndt_subj, zr = zr_subj)
  ) |>
  mutate(response = ifelse(response == 1, "upper", "lower"))

# The 4par summary statistics need both boundaries in every cell: report the
# minority-boundary counts and stop if any cell is too sparse to use.
minority_bias <- data_bias |>
  count(id, stimulus, bias, response) |>
  pivot_wider(names_from = response, values_from = n, values_fill = 0) |>
  mutate(minority = pmin(upper, lower))
cat("Bias dataset -- minority-boundary responses per cell: min",
    min(minority_bias$minority), "| median", median(minority_bias$minority), "\n")
stopifnot(min(minority_bias$minority) >= 15)

write_csv(data_bias, here("data", "simulated_rt_data_bias.csv"))
write_csv(bias_params, here("data", "true_parameters_bias.csv"))
write_csv(
  subject_bias |>
    select(id, stimulus, bias, drift_subj, bound_subj, ndt_subj, zr_subj),
  here("data", "true_parameters_bias_subject.csv")
)
