# Tutorial 1: The EZ-Diffusion Model (3par) on accuracy-coded data
# Tests selective influence in the Wagenmakers et al. (2008) lexical decision
# data: stimulus category should affect drift, emphasis condition the bound.

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

# Default priors for this formula, data, and model
default_prior(formula_acc, data_acc_ezdm, model_acc)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

fit_acc_ezdm <- bmm(
  formula = formula_acc,
  data    = data_acc_ezdm,
  model   = model_acc,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
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
summary(fit_acc_ezdm)
rhat_max <- max(brms::rhat(fit_acc_ezdm), na.rm = TRUE)
cat("Max Rhat:", rhat_max, "\n")

## 3.2) Posterior predictive checks --------------------------------------------
# One panel per observable of the likelihood, each showing its distribution
# across design cells: observed (dark) vs. 50 posterior predictive draws (light)
pp_check_vars(fit_acc_ezdm)
set.seed(1)
pp_draws <- sample(ndraws(fit_acc_ezdm), 50)

p_pp_mean <- pp_check(fit_acc_ezdm, resp_var = "mean_rt", draw_ids = pp_draws) +
  labs(title = "Mean RT (s)", subtitle = NULL, x = NULL) + clean_plot()
p_pp_var <- pp_check(fit_acc_ezdm, resp_var = "var_rt", draw_ids = pp_draws) +
  labs(title = "RT variance (s²)", subtitle = NULL, x = NULL) + clean_plot()
# With accuracy coding the upper boundary is "correct"; legend shown once below
p_pp_pc <- pp_check(fit_acc_ezdm, resp_var = "mean_pc", draw_ids = pp_draws) +
  labs(title = "Proportion correct", subtitle = NULL, x = NULL) + clean_plot() +
  guides(colour = "none", fill = "none")

p_pp <- wrap_plots(p_pp_mean, p_pp_var, p_pp_pc, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(here("figures", "11_pp_check.pdf"), p_pp,
       width = 6.5, height = 3.2, device = "pdf")

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

# One row per posterior draw, stacked so both categories share one panel
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

# Bound draws by condition (log link -> exponentiate)
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
# Savage-Dickey equality tests: Evid.Ratio is the BF for equality, so
# 1 / Evid.Ratio is the evidence against it

# Does drift differ between words and nonwords?
hyp_drift <- hypothesis(fit_acc_ezdm,
  "drift_stim_catword = drift_stim_catnonword")
print(hyp_drift)
cat("\nBF against equal drift (word vs. nonword):",
    round(1 / hyp_drift$hypothesis$Evid.Ratio, 1), "\n")

# Does bound differ between accuracy and speed emphasis?
hyp_bound <- hypothesis(fit_acc_ezdm,
  "bound_conditionaccuracy = bound_conditionspeed")
print(hyp_bound)
cat("\nBF against equal bound (accuracy vs. speed):",
    round(1 / hyp_bound$hypothesis$Evid.Ratio, 1), "\n")
