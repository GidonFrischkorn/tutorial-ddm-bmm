#' Prepare Simulated Data for EAM Tutorials
#'
#' Generates trial-level data from known DDM parameters using rddm().
#' Produces two versions:
#'   1) Clean data -- ground truth for parameter recovery
#'   2) Contaminated data -- with injected contaminants from three processes
#'
#' The three contaminant processes are the ones evaluated in the Monte Carlo
#' study of scripts/02_simulation_contaminant_methods.R, so the demonstration
#' and the simulation describe the same phenomena. Each contaminated trial
#' carries a contaminant_type label, which lets Tutorial 0 show where each
#' process falls relative to the core RT distribution.
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
# 3) Inject Contaminants -------------------------------------------------------
###############################################################################!

# Create a contaminated version of the data by replacing 3% of trials with
# each of three contaminant processes (9% in total). The three processes are
# the ones evaluated in the Monte Carlo study (Appendix A) and differ in where
# they fall relative to the core RT distribution:
#
#   a) Leading-edge anticipations: each subject draws a depth x_i ~ U(0.1, 0.5)
#      (shared across conditions); RT ~ U((1 - x_i) ndt_subj, (1 + x_i) ndt_subj)
#      with accuracy at chance -- premature responses that straddle the leading
#      edge of the core RT distribution, so neither hard cutoffs nor the mixture
#      model gets a free case (Lerche, Voss & Nagler, 2017, adapted to centre
#      on ndt; the v7 construct of the Monte Carlo study).
#   b) Delayed start-ups: Uniform(0, 2 s) added to the RT the decision process
#      produced, leaving the response untouched -- an attentional lapse that
#      delays the start of processing, which places part of the contaminants
#      inside the core distribution (Ratcliff & Tuerlinckx, 2002).
#   c) Informationless responses: the trial is redrawn from the same DDM with
#      drift = 0 -- a response given without stimulus information, at ordinary
#      speed and with accuracy at chance, which overlaps the core completely
#      (Vandekerckhove, Tuerlinckx, & Lee, 2008).

contam_rate <- 0.03  # per process, per participant x condition cell

data_contaminated <- data_clean |>
  group_by(id, condition) |>
  mutate(
    n = n(),
    contam_draw = sample(n()),
    k = round(n * contam_rate),
    contaminant_type = case_when(
      contam_draw <= k     ~ "leading-edge anticipation",
      contam_draw <= 2 * k ~ "delayed start-up",
      contam_draw <= 3 * k ~ "informationless",
      .default = "none"
    )
  ) |>
  ungroup()

# Leading-edge anticipations: each subject draws an anticipation depth
# x_i ~ U(0.1, 0.5), shared across conditions; contaminated RTs fall in a
# uniform band around the subject's ndt, response at chance
subj_x <- tibble(id = 1:n_subjects, x_i = runif(n_subjects, 0.1, 0.5))
data_contaminated <- data_contaminated |> left_join(subj_x, by = "id")
is_fast <- data_contaminated$contaminant_type == "leading-edge anticipation"
data_contaminated$rt[is_fast] <- runif(
  sum(is_fast),
  (1 - data_contaminated$x_i[is_fast]) * data_contaminated$ndt_true[is_fast],
  (1 + data_contaminated$x_i[is_fast]) * data_contaminated$ndt_true[is_fast]
)
data_contaminated$response[is_fast] <- rbinom(sum(is_fast), 1, 0.5)

# Delayed start-ups: add a uniform delay, keep the response
is_delay <- data_contaminated$contaminant_type == "delayed start-up"
data_contaminated$rt[is_delay] <- data_contaminated$rt[is_delay] +
  runif(sum(is_delay), 0, 2)

# Informationless responses: redraw the trial with drift = 0
is_guess <- data_contaminated$contaminant_type == "informationless"
redraw <- rddm(sum(is_guess), drift = 0,
               bound = data_contaminated$bound_true[is_guess],
               ndt = data_contaminated$ndt_true[is_guess], zr = 0.5)
data_contaminated$rt[is_guess]       <- redraw$rt
data_contaminated$response[is_guess] <- redraw$response

data_contaminated <- data_contaminated |>
  mutate(is_contaminant = contaminant_type != "none") |>
  select(id, condition, trial, rt, response, is_contaminant, contaminant_type)

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
