#' Prepare Simulated Data for EAM Tutorials
#'
#' Generates trial-level data from known DDM parameters using rddm().
#' Produces two versions:
#'   1) Clean data -- ground truth for parameter recovery
#'   2) Contaminated data -- with injected fast guesses and slow outliers
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

contam_rate_fast <- 0.04
contam_rate_slow <- 0.03

# Compute per-subject RT statistics for slow outlier range
rt_stats <- data_clean |>
  group_by(id) |>
  summarize(
    mean_rt_subj = mean(rt),
    sd_rt_subj   = sd(rt),
    .groups = "drop"
  )

data_contaminated <- data_clean |>
  left_join(rt_stats, by = "id") |>
  group_by(id, condition) |>
  mutate(
    n = n(),
    # Randomly assign contamination status
    contam_draw = sample(n()),
    is_fast_contam = contam_draw <= round(n * contam_rate_fast),
    is_slow_contam = !is_fast_contam &
      contam_draw > (n - round(n * contam_rate_slow))
  ) |>
  ungroup() |>
  mutate(
    # Inject fast guesses: RT below ndt, response at chance
    rt = if_else(is_fast_contam,
                 runif(n(), 0.1, ndt_true),
                 rt),
    response = if_else(is_fast_contam,
                       as.numeric(rbinom(n(), 1, 0.5)),
                       response),
    # Inject slow outliers: RT far above mean
    rt = if_else(is_slow_contam,
                 runif(n(), mean_rt_subj + 2 * sd_rt_subj,
                       mean_rt_subj + 5 * sd_rt_subj),
                 rt),
    is_contaminant = is_fast_contam | is_slow_contam
  ) |>
  select(id, condition, trial, rt, response, is_contaminant)

###############################################################################!
# 4) Compute EZDM Summary Statistics -------------------------------------------
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
# 5) Save Data -----------------------------------------------------------------
###############################################################################!

# Clean trial-level data (ground truth for parameter recovery)
write_csv(
  data_clean |> select(id, condition, trial, rt, response),
  here("data", "simulated_rt_data_clean.csv")
)

# Contaminated trial-level data (for Tutorial 0 bias demonstration)
write_csv(
  data_contaminated,
  here("data", "simulated_rt_data_contaminated.csv")
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
