#' Tutorial 2: Estimating Response Bias with the EZ-Diffusion Model (4-Parameter)
#'
#' Fits the 4-parameter EZ-Diffusion Model to response-coded Wagenmakers et al. (2008)
#' lexical decision data. Upper boundary = "word" response, lower = "nonword".
#' Drift now reflects evidence toward saying "word"; zr captures response bias.
#' Demonstrates the trial-count constraint on boundary-specific statistics.

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
data_resp_ezdm <- read_csv(here("data", "lexical_response_ezdm4par.csv"),
                            show_col_types = FALSE)
data_response <- read_csv(here("data", "lexical_response_triallevel.csv"),
                           show_col_types = FALSE)

# -- Examine boundary trial counts ---------------------------------------------
# The 4par EZ-Diffusion Model needs stable RT moments on BOTH boundaries in every cell;
# check the minority-boundary counts before fitting
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
# Rhat < 1.01 for all parameters indicates converged chains
summary(fit_resp_ezdm)
rhat_max <- max(brms::rhat(fit_resp_ezdm), na.rm = TRUE)
cat("Max Rhat:", rhat_max, "\n")

## 3.2) Posterior predictive checks --------------------------------------------
# One check per observable of the 4par likelihood: mean RT and RT variance at
# each boundary, and the proportion of upper-boundary ("word") responses. Each
# panel compares how the statistic is distributed across design cells: observed
# (dark) against the posterior predictive draws (light), with the same 50 draws
# in every panel. Where a simulated cell has fewer than two responses at a
# boundary, pp_check() drops that posterior draw (with a warning), so sparse
# minority-boundary cells stay in the check. Cells whose observed boundary has
# no usable RT summaries carry a placeholder and are left out of the RT panels
# (bmm #430); they still enter the proportion panel through their counts.
pp_vars <- pp_check_vars(fit_resp_ezdm)
# Short panel titles; with response coding the upper boundary is "word"
pp_labels <- c(mean_rt_upper = "Mean RT (s), word",
               mean_rt_lower = "Mean RT (s), nonword",
               var_rt_upper  = "RT variance (s\u00b2), word",
               var_rt_lower  = "RT variance (s\u00b2), nonword",
               mean_pc       = "Proportion of word responses")
set.seed(1)
pp_draws <- sample(ndraws(fit_resp_ezdm), 50)
p_pp <- map2(pp_vars$resp_var, pp_labels[pp_vars$resp_var], \(v, lbl)
  pp_check(fit_resp_ezdm, resp_var = v, draw_ids = pp_draws) +
    labs(title = lbl, subtitle = NULL, x = NULL) + clean_plot() +
    # short titles at 9 pt: two panels share the width, and longer titles were
    # clipped at the panel edge in the typeset figure
    theme(plot.title = element_text(size = 9)) +
    # the binned-bars panel carries its own y/y_rep guide; dropping it stops
    # plot_layout(guides = "collect") from printing the legend twice
    (if (v == "mean_pc") guides(colour = "none", fill = "none") else NULL)
) |>
  wrap_plots(ncol = 2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(here("figures", "12_pp_check.pdf"), p_pp,
       width = 6.5, height = 7, device = "pdf")

## 3.3) Parameter estimates ----------------------------------------------------
fixef_est <- fixef(fit_resp_ezdm)
print(round(fixef_est, 3))

cat("\nDrift estimates (response-coded):\n")
cat("  Word:   ", round(fixef_est["drift_stim_catword", "Estimate"], 3), "\n")
cat("  Nonword:", round(fixef_est["drift_stim_catnonword", "Estimate"], 3), "\n")

cat("\nBound estimates:\n")
cat("  Speed:   ", round(exp(fixef_est["bound_conditionspeed", "Estimate"]), 3), "\n")
cat("  Accuracy:", round(exp(fixef_est["bound_conditionaccuracy", "Estimate"]), 3), "\n")

cat("\nRelative starting point (zr):\n")
cat("  Speed:   ", round(plogis(fixef_est["zr_conditionspeed", "Estimate"]), 3), "\n")
cat("  Accuracy:", round(plogis(fixef_est["zr_conditionaccuracy", "Estimate"]), 3), "\n")

## 3.4) Visualize response-coded parameters ------------------------------------

# spread_draws() + pivot_longer(): one row per posterior draw, stacked by
# coefficient (see Tutorial 1 for details)
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

# Savage-Dickey equality tests (see Tutorial 1); the BF against equality is
# 1 / Evid.Ratio

# Test: does drift differ from zero for each stimulus category?
hyp_drift_word <- hypothesis(fit_resp_ezdm, "drift_stim_catword = 0")
print(hyp_drift_word)
cat("\nBF against zero drift (word):",
    round(1 / hyp_drift_word$hypothesis$Evid.Ratio, 1), "\n")

hyp_drift_nw <- hypothesis(fit_resp_ezdm, "drift_stim_catnonword = 0")
print(hyp_drift_nw)
cat("\nBF against zero drift (nonword):",
    round(1 / hyp_drift_nw$hypothesis$Evid.Ratio, 1), "\n")

# Test: does bound differ between accuracy and speed emphasis?
hyp_bound <- hypothesis(fit_resp_ezdm,
  "bound_conditionaccuracy = bound_conditionspeed")
print(hyp_bound)
cat("\nBF against equal bound (accuracy vs. speed):",
    round(1 / hyp_bound$hypothesis$Evid.Ratio, 1), "\n")

# Test: is the starting point biased? zr uses a logit link, so a coefficient
# of 0 corresponds to zr = 0.5 (unbiased)
hyp_zr_acc <- hypothesis(fit_resp_ezdm, "zr_conditionaccuracy = 0")
print(hyp_zr_acc)
cat("\nBF against unbiased starting point (accuracy):",
    round(1 / hyp_zr_acc$hypothesis$Evid.Ratio, 1), "\n")

hyp_zr_speed <- hypothesis(fit_resp_ezdm, "zr_conditionspeed = 0")
print(hyp_zr_speed)
cat("\nBF against unbiased starting point (speed):",
    round(1 / hyp_zr_speed$hypothesis$Evid.Ratio, 1), "\n")
