# Tutorial 3 (continued): The full DDM at the freq specification
# Fits the Diffusion Decision Model to accuracy-coded (zr fixed at 0.5) and
# response-coded (zr estimated) data and compares the estimates with
# Wagenmakers et al. (2008, Table 2).

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

data_accuracy_clean <- read_csv(
  here("data", "lexical_accuracy_triallevel.csv"),
  show_col_types = FALSE
) |>
  mutate(freq = factor(freq, levels = freq_levels))

data_response_clean <- read_csv(
  here("data", "lexical_response_triallevel.csv"),
  show_col_types = FALSE
) |>
  mutate(freq = factor(freq, levels = freq_levels))

# -- Reference parameters (Wagenmakers et al., 2008, rescaled to s = 1) -------
ref_params <- read_csv(here("data", "wagenmakers2008_reference_params.csv"),
                        show_col_types = FALSE)

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

# One model object serves both codings; only the formula and data change
model_ddm <- ddm(rt = "rt", response = "response")

# Accuracy-coded (upper = correct): all drifts are positive
formula_trial_acc <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

# Response-coded (upper = "word"): drift is positive for words, negative for
# nonwords; zr is the bias toward "word" (0.5 = unbiased)
formula_trial_resp <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + condition + (0 + condition | id)
)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

# Each fit takes many hours on this dataset; results are cached via `file`
fit_ddm_acc <- bmm(
  formula = formula_trial_acc,
  data    = data_accuracy_clean,
  model   = model_ddm,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ddm_lexical_acc")
)

fit_ddm_resp <- bmm(
  formula = formula_trial_resp,
  data    = data_response_clean,
  model   = model_ddm,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ddm_lexical_resp")
)

###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) Convergence diagnostics ------------------------------------------------
summary(fit_ddm_acc)
rhat_acc <- max(brms::rhat(fit_ddm_acc), na.rm = TRUE)
cat("Max Rhat (accuracy-coded):", rhat_acc, "\n")

summary(fit_ddm_resp)
rhat_resp <- max(brms::rhat(fit_ddm_resp), na.rm = TRUE)
cat("Max Rhat (response-coded):", rhat_resp, "\n")

## 3.2) Posterior predictive checks --------------------------------------------
# RT, response, and signed RT from one joint simulation; observed (dark) vs.
# simulated (light)
p_pp_acc <- pp_check(fit_ddm_acc, resp_var = "all", ndraws = 50)
ggsave(here("figures", "14_pp_check_acc.pdf"), p_pp_acc,
       width = 6.5, height = 4, device = "pdf")

p_pp_resp <- pp_check(fit_ddm_resp, resp_var = "all", ndraws = 50)
ggsave(here("figures", "14_pp_check_resp.pdf"), p_pp_resp,
       width = 6.5, height = 4, device = "pdf")

## 3.3) Parameter estimates vs. reference values --------------------------------
fe_acc  <- fixef(fit_ddm_acc)
fe_resp <- fixef(fit_ddm_resp)

# Drift by frequency (identity link)
ref_drift <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  select(level, ref_value = value_s1)

cat("\n=== DDM drift by frequency vs. reference (response-coded) ===\n")
cat(sprintf("%-12s %8s %8s\n", "Level", "DDM", "Ref"))
for (lev in freq_levels) {
  ddm_val <- round(fe_resp[paste0("drift_freq", lev), "Estimate"], 3)
  ref_val <- ref_drift$ref_value[ref_drift$level == lev]
  cat(sprintf("%-12s %8.3f %8.3f\n", lev, ddm_val, ref_val))
}

# Bound by condition (log link -> exponentiate)
cat("\n=== DDM bound by condition vs. reference (accuracy-coded) ===\n")
cat("  Speed:   ", round(exp(fe_acc["bound_conditionspeed", "Estimate"]), 3),
    " (ref: 0.84)\n")
cat("  Accuracy:", round(exp(fe_acc["bound_conditionaccuracy", "Estimate"]), 3),
    " (ref: 1.59)\n")

# ndt (log link) and zr (logit link)
cat("\nNdt:", round(exp(fe_acc["ndt_Intercept", "Estimate"]), 3),
    " (ref: 0.408)\n")
cat("\nRelative starting point (zr, response-coded):\n")
cat("  Speed:   ", round(plogis(fe_resp["zr_conditionspeed", "Estimate"]), 3),
    " (ref: 0.476)\n")
cat("  Accuracy:", round(plogis(fe_resp["zr_conditionaccuracy", "Estimate"]), 3),
    " (ref: 0.497)\n")

## 3.4) Visualize parameter estimates ------------------------------------------
draws_drift <- fit_ddm_acc |>
  spread_draws(b_drift_freqhigh, b_drift_freqlow,
               b_drift_freqvery_low, b_drift_freqnonword) |>
  pivot_longer(cols = starts_with("b_drift"),
               names_to = "freq", values_to = "drift") |>
  mutate(freq = str_remove(freq, "b_drift_freq") |>
           factor(levels = freq_levels))

p_drift <- draws_drift |>
  ggplot(aes(x = drift, y = freq, fill = freq)) +
  stat_halfeye(.width = c(.66, .95), alpha = 0.8) +
  scale_fill_manual(values = eam_palette[1:4]) +
  labs(x = expression("Drift ("*delta*") - accuracy-coded"), y = "") +
  guides(fill = "none") +
  clean_plot()

draws_bound <- fit_ddm_acc |>
  spread_draws(b_bound_conditionspeed, b_bound_conditionaccuracy) |>
  pivot_longer(cols = starts_with("b_bound"),
               names_to = "condition", values_to = "bound") |>
  mutate(
    bound = exp(bound),
    condition = str_remove(condition, "b_bound_condition") |> str_to_title()
  )

p_bound <- draws_bound |>
  ggplot(aes(x = bound, y = condition, fill = condition)) +
  stat_halfeye(.width = c(.66, .95), alpha = 0.8) +
  scale_fill_manual(values = eam_palette[3:4]) +
  labs(x = expression("Bound ("*alpha*")"), y = "") +
  guides(fill = "none") +
  clean_plot()

p_params <- p_drift / p_bound + plot_layout(heights = c(1, 1))

ggsave(here("figures", "14_param_estimates.pdf"),
       p_params, width = 6.5, height = 5, device = "pdf")

## 3.5) Hypothesis testing -----------------------------------------------------
# Savage-Dickey tests as in Tutorial 1: 1 / Evid.Ratio is the BF against equality

# Does drift differ between high- and very-low-frequency words?
hyp_drift <- hypothesis(fit_ddm_acc, "drift_freqhigh = drift_freqvery_low",
                        seed = 2026)
print(hyp_drift)
cat("\nBF against equal drift (high vs. very low frequency):",
    round(1 / hyp_drift$hypothesis$Evid.Ratio, 1), "\n")

# Does bound differ between accuracy and speed emphasis?
hyp_bound <- hypothesis(fit_ddm_acc,
  "bound_conditionaccuracy = bound_conditionspeed", seed = 2026)
print(hyp_bound)
cat("\nBF against equal bound (accuracy vs. speed):",
    round(1 / hyp_bound$hypothesis$Evid.Ratio, 1), "\n")
