# Tutorial 2: The EZ-Diffusion Model (4par) on response-coded data
# Upper boundary = "word" response, lower = "nonword": drift is evidence toward
# "word" and zr captures response bias.

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
data_resp_ezdm <- read_csv(here("data", "lexical_response_ezdm4par.csv"),
                            show_col_types = FALSE)
data_response <- read_csv(here("data", "lexical_response_triallevel.csv"),
                           show_col_types = FALSE)

# -- Examine boundary trial counts ---------------------------------------------
# The 4par version needs RT moments at both boundaries in every cell
boundary_counts <- data_response |>
  count(id, condition, stim_cat, response) |>
  pivot_wider(names_from = response, values_from = n, values_fill = 0) |>
  mutate(total = upper + lower, minority = pmin(upper, lower))

cat("Boundary trial counts (response-coded):\n")
boundary_counts |>
  group_by(stim_cat, condition) |>
  summarise(
    median_minority = median(minority),
    min_minority = min(minority),
    .groups = "drop"
  ) |>
  print()

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

model_resp <- ezdm(
  mean_rt  = c("mean_rt_upper", "mean_rt_lower"),
  var_rt   = c("var_rt_upper", "var_rt_lower"),
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "4par"
)

formula_resp_ezdm <- bmf(
  drift ~ 0 + stim_cat + (0 + stim_cat | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + condition + (0 + condition | id)
)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

fit_resp_ezdm <- bmm(
  formula = formula_resp_ezdm,
  data    = data_resp_ezdm,
  model   = model_resp,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_resp")
)

###############################################################################!
# 3) Model Evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) Convergence diagnostics ------------------------------------------------
summary(fit_resp_ezdm)
rhat_max <- max(brms::rhat(fit_resp_ezdm), na.rm = TRUE)
cat("Max Rhat:", rhat_max, "\n")

## 3.2) Posterior predictive checks --------------------------------------------
# Same check as Tutorial 1, now with RT moments at each boundary. Draws with
# fewer than two simulated responses at a boundary are dropped with a warning.
pp_check_vars(fit_resp_ezdm)
set.seed(1)
pp_draws <- sample(ndraws(fit_resp_ezdm), 50)

p_pp_mean_upper <- pp_check(fit_resp_ezdm, resp_var = "mean_rt_upper", draw_ids = pp_draws) +
  labs(title = "Mean RT (s), word", subtitle = NULL, x = NULL) + clean_plot() +
  theme(plot.title = element_text(size = 9))
p_pp_mean_lower <- pp_check(fit_resp_ezdm, resp_var = "mean_rt_lower", draw_ids = pp_draws) +
  labs(title = "Mean RT (s), nonword", subtitle = NULL, x = NULL) + clean_plot() +
  theme(plot.title = element_text(size = 9))
p_pp_var_upper <- pp_check(fit_resp_ezdm, resp_var = "var_rt_upper", draw_ids = pp_draws) +
  labs(title = "RT variance (s²), word", subtitle = NULL, x = NULL) + clean_plot() +
  theme(plot.title = element_text(size = 9))
p_pp_var_lower <- pp_check(fit_resp_ezdm, resp_var = "var_rt_lower", draw_ids = pp_draws) +
  labs(title = "RT variance (s²), nonword", subtitle = NULL, x = NULL) + clean_plot() +
  theme(plot.title = element_text(size = 9))
# With response coding the upper boundary is "word"; legend shown once below
p_pp_pc <- pp_check(fit_resp_ezdm, resp_var = "mean_pc", draw_ids = pp_draws) +
  labs(title = "Proportion of word responses", subtitle = NULL, x = NULL) + clean_plot() +
  theme(plot.title = element_text(size = 9)) +
  guides(colour = "none", fill = "none")

# Row 1: mean RTs and proportion; row 2: RT variances and the legend
p_pp <- wrap_plots(
  p_pp_mean_upper, p_pp_mean_lower, p_pp_pc,
  p_pp_var_upper, p_pp_var_lower, guide_area(),
  design = "ABC\nDEF"
) + plot_layout(guides = "collect")
ggsave(here("figures", "12_pp_check.pdf"), p_pp,
       width = 6.5, height = 4.5, device = "pdf")

## 3.3) Parameter estimates ----------------------------------------------------
fixef_est <- fixef(fit_resp_ezdm)
print(round(fixef_est, 3))

cat("\nDrift estimates (response-coded):\n")
cat("  Word:   ", round(fixef_est["drift_stim_catword", "Estimate"], 3), "\n")
cat("  Nonword:", round(fixef_est["drift_stim_catnonword", "Estimate"], 3), "\n")

cat("\nBound estimates:\n")
cat("  Speed:   ", round(exp(fixef_est["bound_conditionspeed", "Estimate"]), 3), "\n")
cat("  Accuracy:", round(exp(fixef_est["bound_conditionaccuracy", "Estimate"]), 3), "\n")

# zr has a logit link -> plogis
cat("\nRelative starting point (zr):\n")
cat("  Speed:   ", round(plogis(fixef_est["zr_conditionspeed", "Estimate"]), 3), "\n")
cat("  Accuracy:", round(plogis(fixef_est["zr_conditionaccuracy", "Estimate"]), 3), "\n")

## 3.4) Visualize response-coded parameters ------------------------------------
draws_drift <- fit_resp_ezdm |>
  spread_draws(b_drift_stim_catword, b_drift_stim_catnonword) |>
  pivot_longer(cols = starts_with("b_drift"),
               names_to = "stim_cat", values_to = "drift") |>
  mutate(stim_cat = str_remove(stim_cat, "b_drift_stim_cat") |> str_to_title())

p_drift <- draws_drift |>
  ggplot(aes(x = drift, y = stim_cat, fill = stim_cat)) +
  stat_halfeye(.width = c(.66, .95), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = eam_palette[1:2]) +
  labs(x = expression("Drift ("*delta*") - evidence toward 'word'"), y = "") +
  guides(fill = "none") +
  clean_plot()

draws_zr <- fit_resp_ezdm |>
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
  labs(x = expression("Rel. Starting Point ("*zeta*") - bias toward 'word'"), y = "") +
  guides(fill = "none") +
  clean_plot()

p_resp_params <- p_drift / p_zr + plot_layout(heights = c(1, 1))

ggsave(here("figures", "12_resp_params.pdf"),
       p_resp_params, width = 6.5, height = 5, device = "pdf")

## 3.5) Hypothesis testing -----------------------------------------------------
# Savage-Dickey tests as in Tutorial 1: 1 / Evid.Ratio is the BF against equality

# Does drift differ from zero for each stimulus category?
hyp_drift_word <- hypothesis(fit_resp_ezdm, "drift_stim_catword = 0", seed = 2026)
print(hyp_drift_word)
cat("\nBF against zero drift (word):",
    round(1 / hyp_drift_word$hypothesis$Evid.Ratio, 1), "\n")

hyp_drift_nw <- hypothesis(fit_resp_ezdm, "drift_stim_catnonword = 0", seed = 2026)
print(hyp_drift_nw)
cat("\nBF against zero drift (nonword):",
    round(1 / hyp_drift_nw$hypothesis$Evid.Ratio, 1), "\n")

# Does bound differ between accuracy and speed emphasis?
hyp_bound <- hypothesis(fit_resp_ezdm,
  "bound_conditionaccuracy = bound_conditionspeed", seed = 2026)
print(hyp_bound)
cat("\nBF against equal bound (accuracy vs. speed):",
    round(1 / hyp_bound$hypothesis$Evid.Ratio, 1), "\n")

# Is the starting point biased? On the logit scale, 0 means zr = 0.5. Both
# conditions run in one seeded call, as in the manuscript
hyp_zr <- hypothesis(fit_resp_ezdm, c(
  accuracy = "zr_conditionaccuracy = 0",
  speed    = "zr_conditionspeed = 0"
), seed = 2026)
print(hyp_zr)
cat("\nBF against unbiased starting point (accuracy):",
    round(1 / hyp_zr$hypothesis$Evid.Ratio[1], 1), "\n")
cat("BF against unbiased starting point (speed):",
    round(1 / hyp_zr$hypothesis$Evid.Ratio[2], 1), "\n")
