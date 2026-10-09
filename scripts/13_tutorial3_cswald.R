#' Part II (Tutorial 3): The Censored-Shifted Wald at the freq specification
#'
#' Fits both Censored-Shifted Wald versions to Wagenmakers et al. (2008) lexical decision data
#' with the same four-level freq variable on drift as the DDM (script 14), so
#' that Part III can compare the approaches parameter by parameter:
#' - Simple (censored): accuracy-coded, appropriate for low error rates (<20%)
#' - Competing-risks (crisk): response-coded, estimates zr, handles higher errors
#'
#' The composite posterior predictive check of all four trial-level fits is
#' produced by script 15, which loads the DDM fits as well.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ---------------------------------------------------------
chains <- 5
cores  <- 5
# 5 chains x 2 threads = 10 workers. Threading is applied to every model so
# that the fitting times reported in the tutorial stay comparable across the
# three DDM estimation approaches.
threads <- 2
warmup <- 1000
iter   <- 2000
adapt_delta   <- 0.90
max_treedepth <- 12

# -- Load preprocessed data ---------------------------------------------------
# Run scripts/01_prepare_lexical_data.R once first to generate these files
# The CSVs already contain the 4-level freq variable; restore its level order
# after reading (read_csv returns plain character columns).
freq_levels <- c("high", "low", "very_low", "nonword")

data_accuracy_clean <- read_csv(here("data", "lexical_accuracy_triallevel.csv"),
                                 show_col_types = FALSE) |>
  mutate(freq = factor(freq, levels = freq_levels))
data_response_clean <- read_csv(here("data", "lexical_response_triallevel.csv"),
                                 show_col_types = FALSE) |>
  mutate(freq = factor(freq, levels = freq_levels))

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

# Simple version: treats errors as censored correct responses. Use when error
# rates are below ~20% and zr is not needed. At the freq level the
# accuracy-coded error rates run from 1.8% (high-frequency words, accuracy
# emphasis) to 18.7% in seven of eight cells and reach 30.6% for very-low-
# frequency words under speed emphasis; the specification is kept for
# comparability and the manuscript says so. NOTE: the simple version's bound
# is the distance from the starting point to the correct boundary -- half the
# separation ezdm and ddm report -- so it is doubled below.
model_simple <- cswald(
  rt       = "rt",
  response = "response",
  version  = "simple"
)

# Same accuracy-coded formula as the DDM in script 14
formula_trial_acc <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

# Competing-risks (crisk) version: models correct AND error RTs as two racing
# accumulators sharing one ndt. Use at higher error rates, or whenever the
# relative starting point zr is of interest (here: response-coded data).
# On the lexical data (2026-09 refit, vectorized likelihood) its bound agreed
# with the DDM's to within 2%; the ~10% underestimation of the April fits did not recur.
model_crisk <- cswald(
  rt       = "rt",
  response = "response",
  version  = "crisk"
)

# Same response-coded formula as the DDM in script 14
formula_trial_resp <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + condition + (0 + condition | id)
)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

fit_cswald_simple <- bmm(
  formula = formula_trial_acc,
  data    = data_accuracy_clean,
  model   = model_simple,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_cswald_lexical_simple")
)

fit_cswald_crisk <- bmm(
  formula = formula_trial_resp,
  data    = data_response_clean,
  model   = model_crisk,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_cswald_lexical_resp")
)

###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) Convergence diagnostics ------------------------------------------------
# Rhat < 1.01 for all parameters indicates converged chains
summary(fit_cswald_simple)
rhat_simple <- max(brms::rhat(fit_cswald_simple), na.rm = TRUE)
cat("Max Rhat (simple):", rhat_simple, "\n")

summary(fit_cswald_crisk)
rhat_crisk <- max(brms::rhat(fit_cswald_crisk), na.rm = TRUE)
cat("Max Rhat (crisk):", rhat_crisk, "\n")

## 3.2) Posterior predictive checks --------------------------------------------
# All observables of the likelihood (RT, response, signed RT) from one joint
# posterior predictive simulation; observed (dark) vs. simulated (light)
p_pp_simple <- pp_check(fit_cswald_simple, resp_var = "all", ndraws = 50)
ggsave(here("figures", "13_pp_check_simple.pdf"), p_pp_simple,
       width = 6.5, height = 4, device = "pdf")

p_pp_crisk <- pp_check(fit_cswald_crisk, resp_var = "all", ndraws = 50)
ggsave(here("figures", "13_pp_check_crisk.pdf"), p_pp_crisk,
       width = 6.5, height = 4, device = "pdf")

## 3.3) Parameter estimates ----------------------------------------------------
fixef_simple <- fixef(fit_cswald_simple)
print(round(fixef_simple, 3))

# Drift by frequency level. The simple version estimates drift on the LOG
# scale (its single accumulator must drift toward the boundary), so
# back-transform with exp() before reading or comparing the values; the
# crisk version, ezdm, and ddm use an identity link for drift.
cat("\nCensored-Shifted Wald (simple) -- drift estimates by frequency (exp of the log-scale estimate):\n")
for (lev in freq_levels) {
  cat(sprintf("  %-9s %6.3f\n", lev,
              exp(fixef_simple[paste0("drift_freq", lev), "Estimate"])))
}

# Bound by condition (log link -> exponentiate). The simple version estimates
# the distance from the starting point to the correct boundary, which is half
# the boundary separation reported by ezdm and ddm, so multiply by 2 to put it
# on the same scale as the other two approaches.
bound_simple <- exp(fixef_simple[c("bound_conditionspeed",
                                   "bound_conditionaccuracy"), "Estimate"]) * 2
cat("\nCensored-Shifted Wald (simple) -- bound estimates (full-boundary scale):\n")
cat("  Speed:   ", round(bound_simple[1], 3), "\n")
cat("  Accuracy:", round(bound_simple[2], 3), "\n")

fixef_crisk <- fixef(fit_cswald_crisk)

# Relative starting point (logit link -> plogis)
cat("\nCensored-Shifted Wald (crisk) -- relative starting point (zr):\n")
cat("  Speed:   ", round(plogis(fixef_crisk["zr_conditionspeed", "Estimate"]), 3), "\n")
cat("  Accuracy:", round(plogis(fixef_crisk["zr_conditionaccuracy", "Estimate"]), 3), "\n")

## 3.4) Visualize parameter estimates ------------------------------------------

# spread_draws() + pivot_longer(): one row per posterior draw, stacked by
# coefficient (see Tutorial 1 for details)
draws_drift <- fit_cswald_simple |>
  spread_draws(b_drift_freqhigh, b_drift_freqlow,
               b_drift_freqvery_low, b_drift_freqnonword) |>
  pivot_longer(cols = starts_with("b_drift"),
               names_to = "freq", values_to = "drift") |>
  mutate(freq = str_remove(freq, "b_drift_freq") |>
           factor(levels = freq_levels),
         drift = exp(drift))  # log link -> natural scale

p_drift <- draws_drift |>
  ggplot(aes(x = drift, y = freq, fill = freq)) +
  stat_halfeye(.width = c(.66, .95), alpha = 0.8) +
  scale_fill_manual(values = eam_palette[1:4]) +
  labs(x = expression("Drift ("*delta*") - Censored-Shifted Wald (simple)"), y = "") +
  guides(fill = "none") +
  clean_plot()

draws_zr <- fit_cswald_crisk |>
  spread_draws(b_zr_conditionspeed, b_zr_conditionaccuracy) |>
  pivot_longer(cols = starts_with("b_zr"),
               names_to = "condition", values_to = "zr") |>
  mutate(zr = plogis(zr),
         condition = str_remove(condition, "b_zr_condition") |> str_to_title())

p_zr <- draws_zr |>
  ggplot(aes(x = zr, y = condition, fill = condition)) +
  stat_halfeye(.width = c(.66, .95), alpha = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = eam_palette[3:4]) +
  labs(x = expression("Rel. Starting Point ("*zeta*") - Censored-Shifted Wald (crisk)"), y = "") +
  guides(fill = "none") +
  clean_plot()

p_params <- p_drift / p_zr + plot_layout(heights = c(1, 1))

ggsave(here("figures", "13_param_estimates.pdf"),
       p_params, width = 6.5, height = 5, device = "pdf")

## 3.5) Hypothesis testing -----------------------------------------------------

# Savage-Dickey equality tests (see Tutorial 1); the BF against equality is
# 1 / Evid.Ratio. These are the Bayes Factors reported in the manuscript.

# Test: does drift differ between high- and very-low-frequency words?
hyp_drift <- hypothesis(fit_cswald_simple,
  "drift_freqhigh = drift_freqvery_low")
print(hyp_drift)
cat("\nBF against equal drift (high vs. very low frequency):",
    round(1 / hyp_drift$hypothesis$Evid.Ratio, 1), "\n")

# Test: does bound differ between accuracy and speed emphasis?
hyp_bound <- hypothesis(fit_cswald_simple,
  "bound_conditionaccuracy = bound_conditionspeed")
print(hyp_bound)
cat("\nBF against equal bound (accuracy vs. speed):",
    round(1 / hyp_bound$hypothesis$Evid.Ratio, 1), "\n")
