# Simulated data for the parameter recovery studies
# Generates DDM data with known parameters using rddm(): an easy/hard design
# with zr = 0.5 (appendix, script 21) and a stimulus x bias design with
# zr = 0.60/0.40 (Part I, script 22).

###############################################################################!
# 0) R Setup: Packages ---------------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, tidyverse)

set.seed(2025)

###############################################################################!
# 1) Simulation Parameters -----------------------------------------------------
###############################################################################!

# 50 participants, 2 conditions (easy vs. hard), 200 trials per condition
n_subjects    <- 50
n_trials_cond <- 200
conditions    <- c("easy", "hard")

# Population-level parameters: only drift differs between conditions
pop_params <- tibble(
  condition = conditions,
  drift     = c(2.5, 1.0),
  bound     = c(1.5, 1.5),
  ndt       = c(0.3, 0.3),
  zr        = c(0.5, 0.5)
)

# Between-subject variability (SD of random effects)
tau_drift <- 0.5   # SD on identity scale
tau_bound <- 0.15  # SD on log scale
tau_ndt   <- 0.03  # SD on identity scale

###############################################################################!
# 2) Generate Clean Data -------------------------------------------------------
###############################################################################!

# Subject-level parameters: population values plus a random effect
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

# rowwise() simulates each participant x condition with its own parameters;
# reframe() expands each row into n_trials_cond trials
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
# 3) Compute EZ-Diffusion Summary Statistics ----------------------------------
###############################################################################!

# 3par: drift, bound, ndt (zr fixed at 0.5)
data_ezdm_3par <- data_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par", method = "simple"),
    .by = c(id, condition)
  )

# 4par: adds zr via boundary-specific moments
data_ezdm_4par <- data_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par", method = "simple"),
    .by = c(id, condition)
  )

###############################################################################!
# 4) Save Data -----------------------------------------------------------------
###############################################################################!

write_csv(
  data_clean |> select(id, condition, trial, rt, response),
  here("data", "simulated_rt_data_clean.csv")
)

write_csv(data_ezdm_3par, here("data", "simulated_ezdm_3par.csv"))
write_csv(data_ezdm_4par, here("data", "simulated_ezdm_4par.csv"))

# True generating parameters: population and subject level
write_csv(pop_params, here("data", "true_parameters.csv"))
write_csv(
  subject_params |>
    select(id, condition, drift_subj, bound_subj, ndt_subj),
  here("data", "true_parameters_subject.csv")
)

###############################################################################!
# 5) Simulated Data With a Biased Starting Point -------------------------------
###############################################################################!

# Stimulus sets the sign of drift (A: +, B: -); bias sets zr toward the "A"
# (0.60) or "B" (0.40) boundary. Generated after the first dataset on purpose.
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

# Smaller drift spread than above, so every cell keeps enough minority-boundary
# responses for the 4par moments; zr varies on the logit scale
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

# Stop if any cell has too few responses at its minority boundary
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
