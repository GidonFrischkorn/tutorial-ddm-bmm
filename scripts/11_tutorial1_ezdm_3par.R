#' Tutorial 1: Exploring Condition Effects with the EZDM (3-Parameter)
#'
#' Fits the 3-parameter EZDM to accuracy-coded Wagenmakers et al. (2008)
#' lexical decision data. Tests selective influence: stimulus category
#' should affect drift, emphasis condition should affect bound.

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
warmup <- 1000
iter   <- 2000
adapt_delta   <- 0.90
max_treedepth <- 12

# -- Load preprocessed data ---------------------------------------------------
# Run scripts/01_prepare_lexical_data.R once first to generate these files
data_acc_ezdm <- read_csv(here("data", "lexical_accuracy_ezdm3par.csv"),
                           show_col_types = FALSE)

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

model_acc <- ezdm(
  mean_rt  = "mean_rt",
  var_rt   = "var_rt",
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "3par"
)

formula_acc <- bmf(
  drift ~ 0 + stim_cat + (0 + stim_cat | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

fit_acc_ezdm <- bmm(
  formula = formula_acc,
  data    = data_acc_ezdm,
  model   = model_acc,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_acc")
)

###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) Convergence diagnostics ------------------------------------------------
# Rhat < 1.01 for all parameters indicates converged chains
summary(fit_acc_ezdm)
rhat_max <- max(brms::rhat(fit_acc_ezdm), na.rm = TRUE)
cat("Max Rhat:", rhat_max, "\n")

## 3.2) Posterior predictive check ---------------------------------------------
# Overlays data simulated from the posterior (light lines) on the observed
# data (dark line); close overlap indicates adequate model fit
p_pp <- pp_check(fit_acc_ezdm, ndraws = 50) + clean_plot()
ggsave(here("figures", "11_pp_check.pdf"), p_pp,
       width = 6.5, height = 4, device = "pdf")

## 3.3) Parameter estimates ----------------------------------------------------
fixef_est <- fixef(fit_acc_ezdm)
print(round(fixef_est, 3))

# Drift by stimulus category (identity link)
cat("\nDrift estimates:\n")
cat("  Word:   ", round(fixef_est["drift_stim_catword", "Estimate"], 3),
    "[", round(fixef_est["drift_stim_catword", "Q2.5"], 3),
    ",", round(fixef_est["drift_stim_catword", "Q97.5"], 3), "]\n")
cat("  Nonword:", round(fixef_est["drift_stim_catnonword", "Estimate"], 3),
    "[", round(fixef_est["drift_stim_catnonword", "Q2.5"], 3),
    ",", round(fixef_est["drift_stim_catnonword", "Q97.5"], 3), "]\n")

# Bound by condition (log link -> exponentiate)
cat("\nBound estimates:\n")
cat("  Speed:   ", round(exp(fixef_est["bound_conditionspeed", "Estimate"]), 3), "\n")
cat("  Accuracy:", round(exp(fixef_est["bound_conditionaccuracy", "Estimate"]), 3), "\n")

## 3.4) Visualize parameter estimates ------------------------------------------

# spread_draws() (tidybayes) returns one row per posterior draw with the
# requested coefficients as columns; pivot_longer() stacks them so both
# stimulus categories can be plotted in one panel
draws_drift <- fit_acc_ezdm |>
  spread_draws(b_drift_stim_catword, b_drift_stim_catnonword) |>
  pivot_longer(
    cols = starts_with("b_drift"),
    names_to = "stim_cat",
    values_to = "drift"
  ) |>
  mutate(stim_cat = str_remove(stim_cat, "b_drift_stim_cat") |> str_to_title())

p_drift <- draws_drift |>
  ggplot(aes(x = drift, y = stim_cat, fill = stim_cat)) +
  stat_halfeye(.width = c(.66, .95), alpha = 0.8) +
  scale_fill_manual(values = eam_palette[1:2]) +
  labs(x = expression("Drift ("*delta*")"), y = "") +
  guides(fill = "none") +
  clean_plot()

# Extract bound draws by condition (exponentiate)
draws_bound <- fit_acc_ezdm |>
  spread_draws(b_bound_conditionspeed, b_bound_conditionaccuracy) |>
  pivot_longer(
    cols = starts_with("b_bound"),
    names_to = "condition",
    values_to = "bound"
  ) |>
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

ggsave(here("figures", "11_param_estimates.pdf"),
       p_params, width = 6.5, height = 5, device = "pdf")

## 3.5) Hypothesis testing -----------------------------------------------------

# Savage-Dickey equality tests: Evid.Ratio is the BF *for* equality; its
# inverse (reported in the manuscript) quantifies the evidence against equality

# Test: does drift differ between words and nonwords?
hyp_drift <- hypothesis(fit_acc_ezdm,
  "drift_stim_catword = drift_stim_catnonword")
print(hyp_drift)
cat("\nBF against equal drift (word vs. nonword):",
    round(1 / hyp_drift$hypothesis$Evid.Ratio, 1), "\n")

# Test: does bound differ between accuracy and speed emphasis?
hyp_bound <- hypothesis(fit_acc_ezdm,
  "bound_conditionaccuracy = bound_conditionspeed")
print(hyp_bound)
cat("\nBF against equal bound (accuracy vs. speed):",
    round(1 / hyp_bound$hypothesis$Evid.Ratio, 1), "\n")
