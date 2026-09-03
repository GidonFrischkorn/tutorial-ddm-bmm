#' Cross-Model Comparison
#'
#' Fits the EZ-Diffusion comparison models with the 4-level freq specification
#' (matching the DDM from Tutorial 4), then loads all fitted models from
#' Tutorials 1-4 and compares parameter estimates across EZDM, CSWALD, and DDM
#' for both accuracy-coded and response-coded analyses. Overlays reference
#' values from Wagenmakers et al. (2008), Table 2 (converted to s = 1 scaling).

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

# -- Load fitted models from Tutorials 1-4 -------------------------------------
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

# Fit the EZDM with the same 4-level freq variable on drift as the DDM in
# Tutorial 4. This demonstrates that the EZDM produces qualitatively identical
# results at a fraction of the computational cost. Both fits are cached via
# `file`, so re-running this script loads them from output/.

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
# 2) Accuracy-Coded Comparison -------------------------------------------------
###############################################################################!

# EZDM (stim_cat), EZDM freq, CSWALD (stim_cat), DDM (freq)
# Compare drift at both granularities, bound across all models

acc_transforms <- list(bound = exp)

# The cswald "simple" bound is the distance from the starting point to the
# correct boundary, i.e. half the boundary separation that ezdm and ddm
# report. Doubling puts all four models on the full-boundary scale.
cswald_transforms <- list(bound = \(x) exp(x) * 2)

# Bound comparison: all 4 models, on the full-boundary scale
bound_acc <- bind_rows(
  extract_fixef(fit_acc_ezdm, "EZDM", "bound", acc_transforms),
  extract_fixef(fit_acc_ezdm_freq, "EZDM (freq)", "bound", acc_transforms),
  extract_fixef(fit_cswald_simple, "CSWALD simple", "bound", cswald_transforms),
  extract_fixef(fit_ddm_acc, "DDM", "bound", acc_transforms)
) |>
  clean_levels() |>
  mutate(
    model = factor(model, levels = c("EZDM", "EZDM (freq)", "CSWALD simple", "DDM")),
    parameter = "bound"
  )

# Reference values for bound
ref_bound <- ref_params |>
  filter(parameter == "bound") |>
  transmute(level = condition, ref_value = value_s1)

# Drift by stim_cat: EZDM and CSWALD
drift_stim <- bind_rows(
  extract_fixef(fit_acc_ezdm, "EZDM", "drift"),
  extract_fixef(fit_cswald_simple, "CSWALD simple", "drift")
) |> clean_levels()

# Drift by freq: EZDM freq and DDM
drift_freq <- bind_rows(
  extract_fixef(fit_acc_ezdm_freq, "EZDM (freq)", "drift"),
  extract_fixef(fit_ddm_acc, "DDM", "drift")
) |> clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM (freq)", "DDM")))

# Reference values for drift (accuracy-coded: all positive, NW = |NW|)
ref_drift_freq <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(
    level = level,
    ref_value = ifelse(level == "nonword", abs(value_s1), value_s1)
  )

# Plot drift by frequency: EZDM freq vs DDM with reference
p_drift_freq <- ggplot(drift_freq, aes(x = level, y = estimate,
                                        color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  geom_point(data = ref_drift_freq, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 4)]) +
  labs(x = NULL, y = "Estimate",
       title = "Drift by frequency",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot drift by stim_cat: EZDM vs CSWALD
p_drift_stim <- ggplot(drift_stim, aes(x = level, y = estimate,
                                         color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.4), size = 0.5) +
  scale_color_manual(values = eam_palette[c(1, 3)]) +
  labs(x = NULL, y = "Estimate",
       title = "Drift by stim_cat",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot bound comparison with reference
p_bound <- ggplot(bound_acc, aes(x = level, y = estimate,
                                  color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(data = ref_bound, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 2, 3, 4)]) +
  labs(x = NULL, y = "Estimate",
       title = "Bound by condition",
       color = "Model", shape = "Model") +
  annotate("text", x = 0.5, y = Inf, label = "x = Wagenmakers (2008)",
           hjust = 0, vjust = 1.5, size = 2.5) +
  clean_plot()

p_acc <- (p_drift_freq | p_drift_stim) / p_bound +
  plot_layout(heights = c(1, 1))

ggsave(here("figures", "15_cross_model_accuracy.pdf"),
       p_acc, width = 6.5, height = 6)

###############################################################################!
# 3) Response-Coded Comparison -------------------------------------------------
###############################################################################!

resp_transforms <- list(bound = exp, zr = plogis)

# Bound: all 4 models
bound_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm, "EZDM 4par", "bound", resp_transforms),
  extract_fixef(fit_resp_ezdm_freq, "EZDM 4par (freq)", "bound", resp_transforms),
  extract_fixef(fit_cswald_crisk, "CSWALD crisk", "bound", resp_transforms),
  extract_fixef(fit_ddm_resp, "DDM", "bound", resp_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM 4par", "EZDM 4par (freq)",
                                           "CSWALD crisk", "DDM")))

# zr: all 4 models
zr_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm, "EZDM 4par", "zr", resp_transforms),
  extract_fixef(fit_resp_ezdm_freq, "EZDM 4par (freq)", "zr", resp_transforms),
  extract_fixef(fit_cswald_crisk, "CSWALD crisk", "zr", resp_transforms),
  extract_fixef(fit_ddm_resp, "DDM", "zr", resp_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM 4par", "EZDM 4par (freq)",
                                           "CSWALD crisk", "DDM")))

# Reference values for zr
ref_zr <- ref_params |>
  filter(parameter == "zr") |>
  transmute(level = condition, ref_value = value_s1)

# Drift by stim_cat: EZDM 4par and CSWALD crisk
drift_stim_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm, "EZDM 4par", "drift"),
  extract_fixef(fit_cswald_crisk, "CSWALD crisk", "drift")
) |> clean_levels()

# Drift by freq: EZDM 4par freq and DDM
drift_freq_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm_freq, "EZDM 4par (freq)", "drift"),
  extract_fixef(fit_ddm_resp, "DDM", "drift")
) |> clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM 4par (freq)", "DDM")))

ref_drift_resp <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(level = level, ref_value = value_s1)

# Plot response-coded drift by frequency: EZDM freq vs DDM
p_drift_resp <- ggplot(drift_freq_resp, aes(x = level, y = estimate,
                                             color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  geom_point(data = ref_drift_resp, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = eam_palette[c(1, 4)]) +
  labs(x = NULL, y = "Estimate", title = "Drift by frequency",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot response-coded drift by stim_cat: EZDM vs CSWALD
p_drift_stim_resp <- ggplot(
  drift_stim_resp,
  aes(x = level, y = estimate, color = model, shape = model)
) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.4), size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = eam_palette[c(1, 3)]) +
  labs(x = NULL, y = "Estimate",
       title = "Drift by stim_cat",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot bound with reference
p_bound_resp <- ggplot(bound_resp, aes(x = level, y = estimate,
                                        color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(data = ref_bound, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 2, 3, 4)]) +
  labs(x = NULL, y = "Estimate", title = "Bound by condition",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot zr with reference
p_zr <- ggplot(zr_resp, aes(x = level, y = estimate,
                              color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(data = ref_zr, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = eam_palette[c(1, 2, 3, 4)]) +
  labs(x = NULL, y = "Estimate", title = "Rel. starting point by condition",
       color = "Model", shape = "Model") +
  clean_plot()

p_resp <- (p_drift_resp | p_drift_stim_resp) /
  (p_bound_resp | p_zr) +
  plot_layout(heights = c(1, 1))

ggsave(here("figures", "15_cross_model_response.pdf"),
       p_resp, width = 6.5, height = 6)

###############################################################################!
# 4) EZDM vs DDM Figure (Tutorial 4) --------------------------------------------
###############################################################################!

# Three-panel figure for the manuscript's "EZ-Diffusion Model as a Fast
# Alternative" section: drift by frequency, bound by condition, and ndt for
# the accuracy-coded EZDM (freq) vs DDM, with reference values as X markers.

ndt_transforms <- list(ndt = exp)

comp_ndt <- bind_rows(
  extract_fixef(fit_acc_ezdm_freq, "EZDM (freq)", "ndt", ndt_transforms),
  extract_fixef(fit_ddm_acc, "DDM", "ndt", ndt_transforms)
) |>
  mutate(model = factor(model, levels = c("EZDM (freq)", "DDM")),
         level = "ndt")

ref_ndt <- ref_params |>
  filter(parameter == "ndt") |>
  transmute(level = "ndt", ref_value = value_s1)

p_cmp_drift <- p_drift_freq +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cmp_bound <- ggplot(
  bound_acc |> filter(model %in% c("EZDM (freq)", "DDM")),
  aes(x = level, y = estimate, color = model, shape = model)
) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  geom_point(data = ref_bound, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 4)]) +
  labs(x = NULL, y = NULL, title = "Bound by condition",
       color = "Model", shape = "Model") +
  clean_plot()

p_cmp_ndt <- ggplot(comp_ndt, aes(x = level, y = estimate,
                                   color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  geom_point(data = ref_ndt, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 4)]) +
  labs(x = NULL, y = NULL, title = "Non-decision time",
       color = "Model", shape = "Model") +
  clean_plot()

p_ezdm_vs_ddm <- (p_cmp_drift | p_cmp_bound | p_cmp_ndt) +
  plot_layout(guides = "collect", widths = c(2, 1.2, 0.8)) &
  theme(legend.position = "bottom")

ggsave(here("figures", "15_ezdm_vs_ddm_freq.pdf"),
       p_ezdm_vs_ddm, width = 6.5, height = 3.8)

###############################################################################!
# 5) Computational Cost Summary ------------------------------------------------
###############################################################################!

# Wall-clock times read off the archived fits with rstan::get_elapsed_time()
# (slowest chain, since the 5 chains run in parallel). Keep these in step with
# the Computational Cost section of the manuscript.
cat("\n=== Computational Cost Summary ===\n")
cat("Wall-clock (5 chains in parallel); depends on hardware:\n")
cat("  EZDM 3par (accuracy):     21 s (74 s for the freq specification)\n")
cat("  EZDM 4par (response):     43 s (144 s for the freq specification)\n")
cat("  CSWALD simple (accuracy): 4.5 hours\n")
cat("  CSWALD crisk (response):  17.9 hours\n")
cat("  DDM accuracy:             20.6 hours\n")
cat("  DDM response:             17.5 hours\n")
