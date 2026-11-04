# Tutorial 3: The Censored-Shifted Wald at the freq specification
# Fits the simple version (accuracy-coded, for low error rates) and the
# competing-risks version (response-coded, estimates zr) with the same freq
# predictor on drift as the DDM in script 14.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ---------------------------------------------------------
# Every fit uses 2 threads per chain so that fitting times are comparable
chains <- 5
cores  <- 5
threads <- 2
warmup <- 1000
iter   <- 2000
adapt_delta   <- 0.90
max_treedepth <- 12

# -- Load preprocessed data (created by scripts/01_prepare_lexical_data.R) ----
# read_csv() returns characters, so restore the order of the freq levels
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

# Simple version: errors enter as censored correct responses (error rates < ~20%).
# Its bound is the distance from start to the correct boundary, half the DDM's.
model_simple <- cswald(
  rt       = "rt",
  response = "response",
  version  = "simple"
)

formula_trial_acc <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

# Competing-risks version: correct and error responses race to their boundaries,
# for higher error rates or when zr is of interest
model_crisk <- cswald(
  rt       = "rt",
  response = "response",
  version  = "crisk"
)

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
summary(fit_cswald_simple)
rhat_simple <- max(brms::rhat(fit_cswald_simple), na.rm = TRUE)
cat("Max Rhat (simple):", rhat_simple, "\n")

summary(fit_cswald_crisk)
rhat_crisk <- max(brms::rhat(fit_cswald_crisk), na.rm = TRUE)
cat("Max Rhat (crisk):", rhat_crisk, "\n")

## 3.2) Posterior predictive checks --------------------------------------------
# RT, response, and signed RT from one joint simulation; observed (dark) vs.
# simulated (light)
p_pp_simple <- pp_check(fit_cswald_simple, resp_var = "all", ndraws = 50)
ggsave(here("figures", "13_pp_check_simple.pdf"), p_pp_simple,
       width = 6.5, height = 4, device = "pdf")

p_pp_crisk <- pp_check(fit_cswald_crisk, resp_var = "all", ndraws = 50)
ggsave(here("figures", "13_pp_check_crisk.pdf"), p_pp_crisk,
       width = 6.5, height = 4, device = "pdf")

## 3.3) Parameter estimates ----------------------------------------------------
fixef_simple <- fixef(fit_cswald_simple)
print(round(fixef_simple, 3))

# The simple version estimates drift on the log scale -> exp()
cat("\nCensored-Shifted Wald (simple) -- drift estimates by frequency (exp of the log-scale estimate):\n")
for (lev in freq_levels) {
  cat(sprintf("  %-9s %6.3f\n", lev,
              exp(fixef_simple[paste0("drift_freq", lev), "Estimate"])))
}

# Bound: log link -> exp(), then x 2 for the full boundary separation
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
  scale_y_discrete(labels = c(very_low = "very low")) +
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
  labs(x = expression("Relative starting point ("*zeta[r]*"), Censored-Shifted Wald (crisk)"), y = "") +
  guides(fill = "none") +
  clean_plot()

p_params <- p_drift / p_zr + plot_layout(heights = c(1, 1))

ggsave(here("figures", "13_param_estimates.pdf"),
       p_params, width = 6.5, height = 5, device = "pdf")

## 3.5) Hypothesis testing -----------------------------------------------------
# Savage-Dickey tests as in Tutorial 1: 1 / Evid.Ratio is the BF against equality

# Does drift differ between high- and very-low-frequency words?
hyp_drift <- hypothesis(fit_cswald_simple,
  "drift_freqhigh = drift_freqvery_low", seed = 2026)
print(hyp_drift)
cat("\nBF against equal drift (high vs. very low frequency):",
    round(1 / hyp_drift$hypothesis$Evid.Ratio, 1), "\n")

# Does bound differ between accuracy and speed emphasis?
hyp_bound <- hypothesis(fit_cswald_simple,
  "bound_conditionaccuracy = bound_conditionspeed", seed = 2026)
print(hyp_bound)
cat("\nBF against equal bound (accuracy vs. speed):",
    round(1 / hyp_bound$hypothesis$Evid.Ratio, 1), "\n")
