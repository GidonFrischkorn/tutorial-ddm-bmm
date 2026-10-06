#' Cross-Model Comparison
#'
#' Fits the EZ-Diffusion comparison models with the 4-level freq specification
#' (used in Part I's frequency refit and in Part III), builds the composite
#' posterior predictive check for Part II, then compares parameter estimates
#' across EZDM, CSWALD, and DDM -- all at the freq specification -- for both
#' codings. Overlays reference values from Wagenmakers et al. (2008), Table 2
#' (converted to s = 1 scaling).

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))
source(here("functions", "extract_fixef.R"))

# -- Fitting settings (for the EZDM comparison fits below) ---------------------
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

# -- Load preprocessed trial-level data ----------------------------------------
# Run scripts/01_prepare_lexical_data.R once first to generate these files
freq_levels <- c("high", "low", "very_low", "nonword")

data_accuracy_clean <- read_csv(
  here("data", "lexical_accuracy_triallevel.csv"), show_col_types = FALSE
) |>
  mutate(freq = factor(freq, levels = freq_levels))

data_response_clean <- read_csv(
  here("data", "lexical_response_triallevel.csv"), show_col_types = FALSE
) |>
  mutate(freq = factor(freq, levels = freq_levels))

# -- Load fitted models from scripts 11-14 -------------------------------------
fit_acc_ezdm      <- readRDS(here("output", "fit_ezdm_lexical_acc.rds"))
fit_cswald_simple <- readRDS(here("output", "fit_cswald_lexical_simple.rds"))
fit_ddm_acc       <- readRDS(here("output", "fit_ddm_lexical_acc.rds"))

fit_resp_ezdm     <- readRDS(here("output", "fit_ezdm_lexical_resp.rds"))
fit_cswald_crisk  <- readRDS(here("output", "fit_cswald_lexical_resp.rds"))
fit_ddm_resp      <- readRDS(here("output", "fit_ddm_lexical_resp.rds"))

# -- Load reference parameters (Wagenmakers et al., 2008, s = 1) ---------------
ref_params <- read_csv(
  here("data", "wagenmakers2008_reference_params.csv"),
  show_col_types = FALSE
)

###############################################################################!
# 1) EZ-Diffusion Comparison Fits (freq specification) --------------------------
###############################################################################!

# Fit the EZDM with the same 4-level freq variable on drift as the trial-level
# models in scripts 13 and 14. The accuracy-coded fit is the frequency refit
# shown at the end of Tutorial 1. Both fits are cached via `file`, so
# re-running this script loads them from output/.

## 1.1) EZDM 3par with freq (accuracy-coded) ------------------------------------

# maxit = 500 matches script 01. The freq split is finer than the stim_cat
# split used there (130 cells instead of 68), so the EM's iteration cap bites
# harder here, not less.
data_acc_ezdm_freq <- data_accuracy_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par",
                       method = "mixture", distribution = "exgaussian",
                       maxit = 500),
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

fit_acc_ezdm_freq <- bmm(
  formula = formula_ezdm_freq,
  data    = data_acc_ezdm_freq,
  model   = model_ezdm,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_acc_freq")
)

## 1.2) EZDM 4par with freq (response-coded) ------------------------------------

# Stress test: the 4par EZDM requires sufficient trials on *both* boundaries
# in every cell (id x condition x freq). Use the robust method (MAD-based
# estimators) with min_trials = 3 to retain cells with very few minority-
# boundary responses; the mixture method's EM step fails in sparse cells.
data_resp_ezdm_freq <- data_response_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par",
                       method = "robust",
                       min_trials = 3),
    .by = c(id, condition, freq)
  )

# Check trial counts on the minority boundary
minority_counts <- data_resp_ezdm_freq |>
  mutate(n_lower = n_trials - n_upper,
         n_minority = pmin(n_upper, n_lower)) |>
  group_by(condition, freq) |>
  summarise(
    median_minority = median(n_minority),
    min_minority = min(n_minority),
    .groups = "drop"
  )
cat("\n=== EZDM 4par: Trial counts on minority boundary ===\n")
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

fit_resp_ezdm_freq <- bmm(
  formula = formula_ezdm_resp_freq,
  data    = data_resp_ezdm_freq,
  model   = model_ezdm_4par,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_resp_freq")
)

###############################################################################!
# 2) Composite Posterior Predictive Check (Part II figure) ---------------------
###############################################################################!

# One row per trial-level fit, one column per observable of the likelihood:
# the RT distribution, the response proportions, and the signed RT (lower-
# boundary responses negative), which checks choices and RTs jointly. Each row
# uses the same 50 posterior draws in all three panels.
trial_fits <- list(
  "Censored-Shifted Wald (simple), accuracy-coded" = fit_cswald_simple,
  "Full DDM, accuracy-coded"                       = fit_ddm_acc,
  "Censored-Shifted Wald (crisk), response-coded"  = fit_cswald_crisk,
  "Full DDM, response-coded"                       = fit_ddm_resp
)
trial_vars <- c(rt = "Response time", response = "Response (0 = lower, 1 = upper)",
                signed_rt = "Signed response time")
set.seed(1)
p_pp_trial <- imap(trial_fits, \(fit, model_lbl) {
  draws <- sample(ndraws(fit), 50)
  panels <- imap(trial_vars, \(var_lbl, v)
    pp_check(fit, resp_var = v, draw_ids = draws) +
      labs(title = var_lbl, subtitle = NULL) + clean_plot() +
      theme(legend.position = "none", plot.title = element_text(size = 8)) +
      # observed RTs lie below 3 s (the original study's upper bound); the
      # simulated tails would otherwise stretch the axis and flatten the density
      coord_cartesian(xlim = switch(v, rt = c(0, 3), signed_rt = c(-3, 3), NULL))
  )
  wrap_elements(full = grid::textGrob(model_lbl, x = 0, hjust = 0,
                                      gp = grid::gpar(fontsize = 9, fontface = "bold"))) /
    wrap_plots(panels, nrow = 1) +
    plot_layout(heights = c(0.12, 1))
}) |>
  wrap_plots(ncol = 1)
ggsave(here("figures", "15_pp_check_trial_level.pdf"),
       p_pp_trial, width = 6.5, height = 8.5)

###############################################################################!
# 3) Accuracy-Coded Comparison -------------------------------------------------
###############################################################################!

# All three approaches at the freq specification: EZDM (freq), CSWALD simple
# (freq), DDM (freq). Drift by frequency, bound by condition, ndt.

acc_transforms <- list(bound = exp, ndt = exp)

# The cswald "simple" bound is the distance from the starting point to the
# correct boundary, i.e. half the boundary separation that ezdm and ddm
# report. Doubling puts all three models on the full-boundary scale. Its drift
# has a log link (ezdm and ddm: identity), so it is exponentiated as well.
cswald_transforms <- list(drift = exp, bound = \(x) exp(x) * 2, ndt = exp)

acc_levels <- c("EZ-Diffusion", "Censored-Shifted Wald (simple)", "Full DDM")
acc_all <- bind_rows(
  extract_fixef(fit_acc_ezdm_freq, acc_levels[1], c("drift", "bound", "ndt"), acc_transforms),
  extract_fixef(fit_cswald_simple, acc_levels[2], c("drift", "bound", "ndt"), cswald_transforms),
  extract_fixef(fit_ddm_acc,       acc_levels[3], c("drift", "bound", "ndt"), acc_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = acc_levels),
         level = if_else(level == "", "ndt", level))

# Reference values (Wagenmakers et al., 2008; accuracy-coded drift = |drift|)
ref_drift_freq <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(level = level,
            ref_value = ifelse(level == "nonword", abs(value_s1), value_s1))
ref_bound <- ref_params |>
  filter(parameter == "bound") |>
  transmute(level = condition, ref_value = value_s1)
ref_ndt <- ref_params |>
  filter(parameter == "ndt") |>
  transmute(level = "ndt", ref_value = value_s1)

cmp_panel <- function(d, ref, title, palette, y_lab = NULL, dodge = 0.6) {
  ggplot(d, aes(x = level, y = estimate, color = model, shape = model)) +
    geom_pointrange(aes(ymin = lower, ymax = upper),
                    position = position_dodge(width = dodge), size = 0.5) +
    geom_point(data = ref, aes(x = level, y = ref_value),
               inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
    scale_color_manual(values = palette) +
    labs(x = NULL, y = y_lab, title = title, color = "Model", shape = "Model") +
    clean_plot()
}

p_acc_drift <- cmp_panel(filter(acc_all, parameter == "drift") |>
                           mutate(level = factor(level, levels = freq_levels)),
                         ref_drift_freq, "Drift by frequency",
                         eam_palette[c(1, 3, 4)], y_lab = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_acc_bound <- cmp_panel(filter(acc_all, parameter == "bound"), ref_bound,
                         "Bound by condition", eam_palette[c(1, 3, 4)])
# X markers (Wagenmakers et al., 2008) are explained in the figure caption
p_acc_ndt   <- cmp_panel(filter(acc_all, parameter == "ndt"), ref_ndt,
                         "ndt", eam_palette[c(1, 3, 4)])

p_acc <- (p_acc_drift | p_acc_bound | p_acc_ndt) +
  plot_layout(guides = "collect", widths = c(2, 1.2, 0.8)) &
  theme(legend.position = "bottom")

ggsave(here("figures", "15_cross_model_accuracy.pdf"),
       p_acc, width = 6.5, height = 3.8)

###############################################################################!
# 4) Response-Coded Comparison -------------------------------------------------
###############################################################################!

# EZDM 4par (freq), CSWALD crisk (freq), DDM (freq): drift by frequency,
# bound by condition, zr by condition.
# fit_ezdm_lexical_resp_freq converged under the cumulants likelihood (2026-09
# refit; 123 of 136 cells enter it), so it stays in the response-coded panel.
resp_transforms <- list(bound = exp, zr = plogis)
resp_levels <- c("EZ-Diffusion (4par)", "Censored-Shifted Wald (crisk)", "Full DDM")
resp_all <- bind_rows(
  extract_fixef(fit_resp_ezdm_freq, resp_levels[1], c("drift", "bound", "zr"), resp_transforms),
  extract_fixef(fit_cswald_crisk,   resp_levels[2], c("drift", "bound", "zr"), resp_transforms),
  extract_fixef(fit_ddm_resp,       resp_levels[3], c("drift", "bound", "zr"), resp_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = resp_levels))

ref_drift_resp <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(level = level, ref_value = value_s1)
ref_zr <- ref_params |>
  filter(parameter == "zr") |>
  transmute(level = condition, ref_value = value_s1)

p_resp_drift <- cmp_panel(filter(resp_all, parameter == "drift") |>
                            mutate(level = factor(level, levels = freq_levels)),
                          ref_drift_resp, "Drift by frequency",
                          eam_palette[c(1, 3, 4)], y_lab = "Estimate") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_resp_bound <- cmp_panel(filter(resp_all, parameter == "bound"), ref_bound,
                          "Bound by condition", eam_palette[c(1, 3, 4)])
p_resp_zr <- cmp_panel(filter(resp_all, parameter == "zr"), ref_zr,
                       "Rel. starting point", eam_palette[c(1, 3, 4)]) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.3)

p_resp <- (p_resp_drift | p_resp_bound | p_resp_zr) +
  plot_layout(guides = "collect", widths = c(2, 1.2, 1)) &
  theme(legend.position = "bottom")

ggsave(here("figures", "15_cross_model_response.pdf"),
       p_resp, width = 6.5, height = 3.8)

###############################################################################!
# 5) Computational Cost Summary ------------------------------------------------
###############################################################################!

# Wall-clock of the slowest chain, read from the cached fits with
# rstan::get_elapsed_time(); the manuscript's setup chunk reads the same
# quantity, so nothing here is typed by hand.
cost <- tibble(
  fit = c("EZDM 3par (stim_cat)", "EZDM 3par (freq)", "EZDM 4par (stim_cat)",
          "EZDM 4par (freq)", "CSWALD simple", "CSWALD crisk",
          "DDM accuracy", "DDM response"),
  hours = purrr::map_dbl(
    list(fit_acc_ezdm, fit_acc_ezdm_freq, fit_resp_ezdm, fit_resp_ezdm_freq,
         fit_cswald_simple, fit_cswald_crisk, fit_ddm_acc, fit_ddm_resp),
    \(f) max(rowSums(rstan::get_elapsed_time(f$fit))) / 3600
  )
)
cat("\n=== Computational Cost Summary (wall-clock, slowest chain) ===\n")
print(cost |> mutate(hours = round(hours, 3)), n = Inf)
