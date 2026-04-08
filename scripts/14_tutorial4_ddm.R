#' Tutorial 4: The DDM as Gold Standard
#'
#' Fits the full Diffusion Decision Model to Wagenmakers et al. (2008) lexical
#' decision data. Uses the 4-level frequency variable on drift to match the
#' original paper's specification. Fits accuracy-coded (all drifts positive,
#' zr fixed at 0.5) and response-coded (drift sign reflects stimulus type,
#' zr estimated) variants.
#'
#' Reference values (s = 1) from Wagenmakers et al. (2008), Table 2:
#'   bound: accuracy = 1.59, speed = 0.84
#'   drift: HF = 4.52, LF = 2.75, VLF = 1.80, NW = -2.64 (response-coded)
#'   ndt:   0.408 (constant)
#'   zr:    accuracy = 0.497, speed = 0.476 (response-coded)

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms, rtdists,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ---------------------------------------------------------
chains <- 5
cores  <- 5
warmup <- 1000
iter   <- 2000
adapt_delta   <- 0.90
max_treedepth <- 12

# -- Load preprocessed data ---------------------------------------------------
# Run scripts/prepare_lexical_data.R once first to generate these files
data_accuracy_clean <- read_csv(
  here("data", "lexical_accuracy_triallevel.csv"),
  show_col_types = FALSE
) |>
  mutate(
    freq = factor(
      ifelse(stim_cat == "nonword", "nonword", as.character(frequency)),
      levels = c("high", "low", "very_low", "nonword")
    )
  )

data_response_clean <- read_csv(
  here("data", "lexical_response_triallevel.csv"),
  show_col_types = FALSE
) |>
  mutate(
    freq = factor(
      ifelse(stim_cat == "nonword", "nonword", as.character(frequency)),
      levels = c("high", "low", "very_low", "nonword")
    )
  )

# -- Load reference parameters for comparison ---------------------------------
ref_params <- read_csv(here("data", "wagenmakers2008_reference_params.csv"),
                        show_col_types = FALSE)

###############################################################################!
# 1) DDM on Accuracy-Coded Data -----------------------------------------------
###############################################################################!

# In accuracy coding (upper = correct), ALL drifts are positive: they measure
# evidence toward the correct response regardless of stimulus type.
# Expected rank order: HF (4.52) > LF (2.75) > NW (2.64) > VLF (1.80)
# Note: NW > VLF because nonwords are easier to correctly reject than
# very-low-frequency words are to correctly accept.

model_ddm <- ddm(rt = "rt", response = "response")

formula_ddm_acc <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

fit_ddm_acc <- bmm(
  formula = formula_ddm_acc,
  data    = data_accuracy_clean,
  model   = model_ddm,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ddm_lexical_acc")
)

## Convergence & estimates -----------------------------------------------------
summary(fit_ddm_acc)
rhat_acc <- max(brms::rhat(fit_ddm_acc), na.rm = TRUE)

## Posterior predictive check ---------------------------------------------------
p_pp_acc <- pp_check(fit_ddm_acc, ndraws = 50) + clean_plot()
ggsave(here("figures", "14_pp_check_acc.pdf"), p_pp_acc,
       width = 6.5, height = 4, device = "pdf")


###############################################################################!
# 2) DDM on Response-Coded Data ------------------------------------------------
###############################################################################!

# In response coding (upper = "word"), drift sign reflects stimulus type:
# positive for words (evidence toward "word"), negative for nonwords.
# Expected: HF = 4.52, LF = 2.75, VLF = 1.80, NW = -2.64
# zr captures bias toward saying "word" (0.5 = unbiased)
# Expected: accuracy = 0.497, speed = 0.476 (very small shift)

formula_ddm_resp <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + condition + (0 + condition | id)
)

fit_ddm_resp <- bmm(
  formula = formula_ddm_resp,
  data    = data_response_clean,
  model   = model_ddm,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ddm_lexical_resp")
)

## Convergence & estimates -----------------------------------------------------
summary(fit_ddm_resp)
rhat_resp <- max(brms::rhat(fit_ddm_resp), na.rm = TRUE)

## Posterior predictive check ---------------------------------------------------
p_pp_resp <- pp_check(fit_ddm_resp, ndraws = 50) + clean_plot()
ggsave(here("figures", "14_pp_check_resp.pdf"), p_pp_resp,
       width = 6.5, height = 4, device = "pdf")

###############################################################################!
# 3) EZDM with Same freq Specification ----------------------------------------
###############################################################################!

# Fit the EZDM 3par with the same 4-level freq variable on drift to compare
# against the DDM. This demonstrates that the EZDM produces qualitatively
# identical results at a fraction of the computational cost.

# -- Compute EZDM summary stats by freq (not stim_cat)
data_acc_ezdm_freq <- data_accuracy_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par",
                       method = "mixture", distribution = "exgaussian"),
    .by = c(id, condition, freq)
  )

model_ezdm <- ezdm(
  mean_rt  = "mean_rt",
  var_rt   = "var_rt",
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "3par"
)

formula_ezdm_freq <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

fit_ezdm_freq <- bmm(
  formula = formula_ezdm_freq,
  data    = data_acc_ezdm_freq,
  model   = model_ezdm,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_acc_freq")
)

## Compare EZDM vs DDM --------------------------------------------------------
cat("\n=== EZDM vs DDM: Drift by frequency (accuracy-coded) ===\n")
fe_ezdm <- fixef(fit_ezdm_freq)
fe_ddm  <- fixef(fit_ddm_acc)

cat(sprintf("%-12s %8s %8s\n", "Level", "EZDM", "DDM"))
for (lev in c("high", "low", "very_low", "nonword")) {
  ezdm_val <- round(fe_ezdm[paste0("drift_freq", lev), "Estimate"], 3)
  ddm_val  <- round(fe_ddm[paste0("drift_freq", lev), "Estimate"], 3)
  cat(sprintf("%-12s %8.3f %8.3f\n", lev, ezdm_val, ddm_val))
}

cat("\n=== EZDM vs DDM: Bound by condition ===\n")
cat(sprintf("%-12s %8s %8s\n", "Condition", "EZDM", "DDM"))
for (cond in c("accuracy", "speed")) {
  ezdm_val <- round(exp(fe_ezdm[paste0("bound_condition", cond), "Estimate"]), 3)
  ddm_val  <- round(exp(fe_ddm[paste0("bound_condition", cond), "Estimate"]), 3)
  cat(sprintf("%-12s %8.3f %8.3f\n", cond, ezdm_val, ddm_val))
}

cat("\n=== EZDM vs DDM: Ndt ===\n")
cat("  EZDM:", round(exp(fe_ezdm["ndt_Intercept", "Estimate"]), 3), "\n")
cat("  DDM: ", round(exp(fe_ddm["ndt_Intercept", "Estimate"]), 3), "\n")

###############################################################################!
# 4) EZDM 4par Response-Coded with freq ----------------------------------------
###############################################################################!

# Fit the EZDM 4par with response coding and the 4-level freq variable.
# This is a stress test: the 4par EZDM requires sufficient trials on *both*
# boundaries in every cell (id x condition x freq). With 4 frequency levels,
# some cells (e.g., nonwords under accuracy emphasis) may have very few
# minority-boundary responses.

# -- Compute EZDM 4par summary stats by freq (response-coded)
# Use robust method (MAD-based estimators) with min_trials = 3 to retain cells
# with very few minority-boundary responses. The mixture method requires enough
# trials to fit an EM model, which fails in sparse cells. The robust method
# provides reliable moment estimates even with small n.
data_resp_ezdm_freq <- data_response_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par",
                       method = "robust",
                       min_trials = 3),
    .by = c(id, condition, freq)
  )

# Check trial counts on minority boundary
cat("\n=== EZDM 4par: Trial counts on minority boundary ===\n")
minority_counts <- data_resp_ezdm_freq |>
  mutate(n_lower = n_trials - n_upper,
         n_minority = pmin(n_upper, n_lower)) |>
  group_by(condition, freq) |>
  summarise(
    median_minority = median(n_minority),
    min_minority = min(n_minority),
    .groups = "drop"
  )
print(minority_counts)

model_ezdm_4par <- ezdm(
  mean_rt  = c("mean_rt_upper", "mean_rt_lower"),
  var_rt   = c("var_rt_upper", "var_rt_lower"),
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "4par"
)

formula_ezdm_resp_freq <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + condition + (0 + condition | id)
)

fit_ezdm_resp_freq <- bmm(
  formula = formula_ezdm_resp_freq,
  data    = data_resp_ezdm_freq,
  model   = model_ezdm_4par,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_resp_freq")
)

## Compare EZDM 4par vs DDM (response-coded) -----------------------------------
cat("\n=== EZDM 4par vs DDM: Drift by frequency (response-coded) ===\n")
fe_ezdm_resp <- fixef(fit_ezdm_resp_freq)
fe_ddm_resp  <- fixef(fit_ddm_resp)

cat(sprintf("%-12s %8s %8s\n", "Level", "EZDM 4par", "DDM"))
for (lev in c("high", "low", "very_low", "nonword")) {
  ezdm_val <- round(fe_ezdm_resp[paste0("drift_freq", lev), "Estimate"], 3)
  ddm_val  <- round(fe_ddm_resp[paste0("drift_freq", lev), "Estimate"], 3)
  cat(sprintf("%-12s %9.3f %8.3f\n", lev, ezdm_val, ddm_val))
}

cat("\n=== EZDM 4par vs DDM: Relative starting point ===\n")
cat(sprintf("%-12s %8s %8s\n", "Condition", "EZDM 4par", "DDM"))
for (cond in c("accuracy", "speed")) {
  ezdm_val <- round(plogis(fe_ezdm_resp[paste0("zr_condition", cond), "Estimate"]), 3)
  ddm_val  <- round(plogis(fe_ddm_resp[paste0("zr_condition", cond), "Estimate"]), 3)
  cat(sprintf("%-12s %9.3f %8.3f\n", cond, ezdm_val, ddm_val))
}
