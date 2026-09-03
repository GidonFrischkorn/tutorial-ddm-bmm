#' Systematic Comparison of Contaminant-Flagging Settings
#'
#' Companion to Tutorial 0 (scripts/10_tutorial0_data_preprocessing.R).
#' Systematically compares the settings of bmm's preprocessing tools on the
#' simulated dataset with known contaminants:
#'   1. Core distribution: exgaussian, lognormal, invgaussian
#'   2. Contaminant bounds: data-driven, fixed wide, fixed narrow
#'   3. Hard pre-filtering: none vs. RT < 150ms / > 5000ms
#' plus threshold-based vs. probabilistic trial exclusion and the
#' validate_fast_guesses() check.
#'
#' This is methods-evaluation material for readers who want to know WHICH
#' settings to use; the tutorial itself only demonstrates the defaults.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, tidyverse, patchwork)
source(here("functions", "clean_plot.R"))

data_contaminated <- read_csv(
  here("data", "simulated_rt_data_contaminated.csv"),
  show_col_types = FALSE
)
data_clean <- read_csv(
  here("data", "simulated_rt_data_clean.csv"),
  show_col_types = FALSE
)

# -- Baseline summary statistics (as in Tutorial 0, section 1.2) ---------------
sumstats_contaminated <- data_contaminated |>
  summarize(
    ezdm_summary_stats(rt, response,
                       version = "3par", method = "simple"),
    .by = c(id, condition)
  ) |>
  mutate(source = "Contaminated")

sumstats_clean <- data_clean |>
  summarize(
    ezdm_summary_stats(rt, response,
                       version = "3par", method = "simple"),
    .by = c(id, condition)
  ) |>
  mutate(source = "Clean")

# -- Factor definitions ---------------------------------------------------------

distributions <- c("exgaussian", "lognormal", "invgaussian")

# Contaminant bound settings to compare:
# - "data_driven": c("min", "max") -- adapts to observed RT range per cell
# - "fixed_wide": c(0.05, 10) -- generous fixed bounds covering any RT
# - "fixed_narrow": c(0.1, 3) -- tighter fixed bounds based on typical RT range
bound_settings <- list(
  data_driven  = c("min", "max"),
  fixed_wide   = c(0.05, 10),
  fixed_narrow = c(0.1, 3)
)

###############################################################################!
# 1) Hard Pre-Filtering: RT < 150ms or RT > 5000ms ------------------------------
###############################################################################!

# A common practice is to remove implausibly fast responses (< 150ms,
# faster than any evidence accumulation process) and extremely slow
# responses (> 5000ms, likely attentional lapses). We evaluate whether
# this pre-filtering step improves contaminant detection.

hard_lower <- 0.15  # 150ms in seconds
hard_upper <- 5.0   # 5000ms in seconds

data_prefiltered <- data_contaminated |>
  filter(rt >= hard_lower & rt <= hard_upper)

cat("\nHard filtering removed",
    nrow(data_contaminated) - nrow(data_prefiltered), "of",
    nrow(data_contaminated), "trials (",
    round(100 * (1 - nrow(data_prefiltered) / nrow(data_contaminated)), 1),
    "%)\n")

# Check how many actual contaminants were caught by hard filtering
hard_filter_performance <- data_contaminated |>
  mutate(removed = rt < hard_lower | rt > hard_upper) |>
  summarize(
    n_removed = sum(removed),
    contam_removed = sum(removed & is_contaminant),
    clean_removed = sum(removed & !is_contaminant),
    .by = condition
  )
cat("\nHard filter catches:\n")
print(hard_filter_performance)

# Define pre-filtering factor
prefilter_settings <- list(
  none = data_contaminated,
  hard = data_prefiltered
)

###############################################################################!
# 2) Contaminant Flagging: Distribution x Bounds x Pre-Filtering ----------------
###############################################################################!

# Full factorial: 3 distributions x 3 bounds x 2 pre-filter = 18 combinations
flag_results <- expand_grid(
  dist = distributions,
  bound_name = names(bound_settings),
  prefilter = names(prefilter_settings)
) |>
  pmap(function(dist, bound_name, prefilter) {
    bounds <- bound_settings[[bound_name]]
    input_data <- prefilter_settings[[prefilter]]
    cat("Flagging:", dist, "x", bound_name, "x prefilter =", prefilter, "\n")
    input_data |>
      group_by(id, condition) |>
      mutate(
        contam_prob = flag_contaminant_rts(
          rt, distribution = dist, contaminant_bound = bounds,
          maxit = 500
        )
      ) |>
      ungroup() |>
      mutate(distribution = dist, bound_setting = bound_name,
             prefilter = prefilter)
  }) |>
  bind_rows()

# Classification performance across all three factors
# For the pre-filtered condition, trials removed by hard filtering count as
# true positives (if contaminant) or false positives (if clean) to give a
# fair total-pipeline comparison
classification_full <- flag_results |>
  mutate(flagged = contam_prob > 0.5) |>
  summarize(
    true_pos  = sum(is_contaminant & flagged, na.rm = TRUE),
    false_pos = sum(!is_contaminant & flagged, na.rm = TRUE),
    true_neg  = sum(!is_contaminant & !flagged, na.rm = TRUE),
    false_neg = sum(is_contaminant & !flagged, na.rm = TRUE),
    n_na      = sum(is.na(flagged)),
    .by = c(distribution, bound_setting, prefilter)
  ) |>
  mutate(
    # For hard-prefiltered: add back the trials removed by hard filtering
    true_pos = if_else(
      prefilter == "hard",
      true_pos + hard_filter_performance$contam_removed[1] +
        hard_filter_performance$contam_removed[2],
      true_pos
    ),
    false_pos = if_else(
      prefilter == "hard",
      false_pos + hard_filter_performance$clean_removed[1] +
        hard_filter_performance$clean_removed[2],
      false_pos
    ),
    sensitivity = true_pos / (true_pos + false_neg),
    specificity = true_neg / (true_neg + false_pos),
    balanced_acc = (sensitivity + specificity) / 2,
    precision = true_pos / (true_pos + false_pos),
    f1 = 2 * precision * sensitivity / (precision + sensitivity)
  )

cat("\n=== Classification: Distribution x Bounds x Pre-filter ===\n\n")
classification_full |>
  select(distribution, bound_setting, prefilter, sensitivity,
         specificity, balanced_acc, f1, n_na) |>
  mutate(across(where(is.numeric) & !n_na, \(x) round(x, 3))) |>
  arrange(desc(f1)) |>
  print(n = Inf)

###############################################################################!
# 3) Validating Flagged Contaminants with validate_fast_guesses() ---------------
###############################################################################!

# validate_fast_guesses() provides a principled check that trials flagged by
# flag_contaminant_rts() actually show random guessing behavior (accuracy at
# chance). It uses a Bayesian beta-binomial model to compute a Bayes Factor
# for H0: p(upper) = 0.5 among flagged fast trials.

# Use the default flagging (exgaussian, data-driven bounds, no pre-filter)
flag_default <- flag_results |>
  filter(distribution == "exgaussian",
         bound_setting == "data_driven",
         prefilter == "none")

validation <- flag_default |>
  group_by(id, condition) |>
  reframe(
    validate_fast_guesses(
      contam_flag = contam_prob > 0.5,
      rt_data = rt,
      response = response
    ) |> as_tibble()
  )

cat("\n=== validate_fast_guesses: Evidence for random guessing ===\n\n")
cat("Distribution of BF evidence categories:\n")
validation |>
  count(bf_evidence) |>
  arrange(desc(n)) |>
  print()

cat("\nMedian BF (H0: guessing):", round(median(validation$bf_01,
                                                 na.rm = TRUE), 2), "\n")
cat("% cells where 0.5 in HDI:",
    round(mean(validation$guess_in_hdi, na.rm = TRUE) * 100, 1), "%\n")
cat("Mean accuracy of flagged fast trials:",
    round(mean(validation$prop_upper, na.rm = TRUE), 3), "\n")

###############################################################################!
# 4) Mixture-Method Aggregation: Distribution x Bounds x Pre-Filter -------------
###############################################################################!

# Full factorial for the mixture method in ezdm_summary_stats()
sumstats_by_settings <- expand_grid(
  dist = distributions,
  bound_name = names(bound_settings),
  prefilter = names(prefilter_settings)
) |>
  pmap(function(dist, bound_name, prefilter) {
    bounds <- bound_settings[[bound_name]]
    input_data <- prefilter_settings[[prefilter]]
    cat("Mixture:", dist, "x", bound_name, "x", prefilter, "\n")
    input_data |>
      summarize(
        ezdm_summary_stats(rt, response,
                           version = "3par", method = "mixture",
                           distribution = dist,
                           contaminant_bound = bounds,
                           maxit = 500),
        .by = c(id, condition)
      ) |>
      mutate(distribution = dist, bound_setting = bound_name,
             prefilter = prefilter, approach = "mixture")
  }) |>
  bind_rows()

###############################################################################!
# 5) Filter-Then-Aggregate: Threshold vs Probabilistic Flagging -----------------
###############################################################################!

# Two approaches to using contaminant probabilities for trial exclusion:
#
# A) Threshold-based: exclude trials with contam_prob > 0.5
#    - Deterministic: same result every time
#    - Binary decision boundary at 0.5
#    - Trials near the boundary are treated identically to clear cases
#
# B) Probabilistic: exclude trials with probability = contam_prob
#    - Stochastic: draws a Bernoulli for each trial using contam_prob
#    - Trials with high contam_prob are very likely excluded,
#      trials with moderate contam_prob have proportional chance
#    - Better reflects uncertainty in the classification
#    - Results vary across runs (set seed for reproducibility)
#
# We compare both approaches against the mixture method (which avoids
# binary exclusion entirely by down-weighting in the aggregation step).

# Threshold-based filter: contam_prob > 0.5
sumstats_threshold <- flag_default |>
  filter(contam_prob <= 0.5) |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par", method = "simple"),
    .by = c(id, condition)
  ) |>
  mutate(approach = "filter_threshold")

# Probabilistic filter: rbinom(1, 1, contam_prob) for each trial
set.seed(2025)
sumstats_probabilistic <- flag_default |>
  mutate(exclude = rbinom(n(), size = 1, prob = contam_prob)) |>
  filter(exclude == 0) |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par", method = "simple"),
    .by = c(id, condition)
  ) |>
  mutate(approach = "filter_probabilistic")

# Mixture method (already computed above for exgaussian, data-driven, no prefilter)
sumstats_mixture_default <- sumstats_by_settings |>
  filter(distribution == "exgaussian",
         bound_setting == "data_driven",
         prefilter == "none") |>
  mutate(approach = "mixture")

# Compare all three approaches against clean ground truth
filter_comparison <- bind_rows(
  sumstats_threshold,
  sumstats_probabilistic,
  sumstats_mixture_default
) |>
  left_join(
    sumstats_clean |>
      select(id, condition, mean_rt, var_rt) |>
      rename(mean_rt_clean = mean_rt, var_rt_clean = var_rt),
    by = c("id", "condition")
  ) |>
  mutate(
    bias_mean_rt = mean_rt - mean_rt_clean,
    bias_var_rt  = var_rt - var_rt_clean
  )

filter_comparison_summary <- filter_comparison |>
  summarize(
    mae_mean_rt = median(abs(bias_mean_rt)),
    mae_var_rt  = median(abs(bias_var_rt)),
    .by = c(approach, condition)
  )

cat("\n=== Flagging Approach Comparison (exGaussian, data-driven bounds) ===\n")
cat("MAE for mean RT and var RT by approach and condition:\n\n")
filter_comparison_summary |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(condition, mae_var_rt) |>
  print(n = Inf)

###############################################################################!
# 6) Full Factorial Filter-Then-Aggregate ---------------------------------------
###############################################################################!

# Use threshold-based filtering across all settings for the comprehensive
# comparison (the manuscript discusses the threshold vs probabilistic
# distinction separately)
sumstats_filtered_by_settings <- expand_grid(
  dist = distributions,
  bound_name = names(bound_settings),
  prefilter = names(prefilter_settings)
) |>
  pmap(function(dist, bound_name, prefilter) {
    flag_results |>
      filter(distribution == dist,
             bound_setting == bound_name,
             prefilter == !!prefilter,
             contam_prob <= 0.5) |>
      summarize(
        ezdm_summary_stats(rt, response,
                           version = "3par", method = "simple"),
        .by = c(id, condition)
      ) |>
      mutate(distribution = dist, bound_setting = bound_name,
             prefilter = prefilter, approach = "filter")
  }) |>
  bind_rows()

# Combine all approaches
sumstats_all_approaches <- bind_rows(
  sumstats_by_settings,
  sumstats_filtered_by_settings
)

# Compute bias relative to clean ground truth for all combinations
sumstats_all_bias <- sumstats_all_approaches |>
  left_join(
    sumstats_clean |>
      select(id, condition, mean_rt, var_rt) |>
      rename(mean_rt_clean = mean_rt, var_rt_clean = var_rt),
    by = c("id", "condition")
  ) |>
  mutate(
    bias_mean_rt = mean_rt - mean_rt_clean,
    bias_var_rt  = var_rt - var_rt_clean
  )

# Summary: median absolute bias by all three factors
bias_summary <- sumstats_all_bias |>
  summarize(
    mae_mean_rt = median(abs(bias_mean_rt), na.rm = TRUE),
    mae_var_rt  = median(abs(bias_var_rt), na.rm = TRUE),
    .by = c(approach, distribution, bound_setting, prefilter, condition)
  )

cat("\n=== Summary Statistic Recovery: MAE (All Factors) ===\n")
cat("Top 10 by var_rt MAE per condition:\n\n")
bias_summary |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(condition, mae_var_rt) |>
  group_by(condition) |>
  slice_head(n = 10) |>
  ungroup() |>
  print(n = Inf)

###############################################################################!
# 7) Comprehensive Comparison ---------------------------------------------------
###############################################################################!

## 7.1) Identify best settings --------------------------------------------------

# Rank all distribution x bound x prefilter combinations by var_rt recovery
best_settings <- bias_summary |>
  arrange(condition, approach, mae_var_rt) |>
  group_by(condition, approach) |>
  slice(1) |>
  ungroup()

cat("\n=== Best Settings per Approach x Condition (by var_rt MAE) ===\n\n")
best_settings |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  print(n = Inf)

# Effect of pre-filtering (averaged across other settings)
prefilter_effect <- bias_summary |>
  summarize(
    mae_var_rt = median(mae_var_rt, na.rm = TRUE),
    .by = c(prefilter, condition)
  )
cat("\n=== Pre-filtering Effect (median MAE var_rt) ===\n\n")
prefilter_effect |>
  pivot_wider(names_from = prefilter, values_from = mae_var_rt) |>
  mutate(improvement = round(1 - hard / none, 3)) |>
  print()

## 7.2) Visualization: classification heatmap (distribution x bounds x prefilter)

p_classification <- classification_full |>
  mutate(
    distribution = factor(distribution,
      levels = c("exgaussian", "lognormal", "invgaussian"),
      labels = c("Ex-Gaussian", "Lognormal", "Inv-Gaussian")),
    bound_setting = factor(bound_setting,
      levels = c("data_driven", "fixed_wide", "fixed_narrow"),
      labels = c("Data-driven", "Fixed wide", "Fixed narrow")),
    prefilter = factor(prefilter,
      levels = c("none", "hard"),
      labels = c("No pre-filter", "Hard pre-filter"))
  ) |>
  ggplot(aes(x = distribution, y = bound_setting, fill = f1)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(f1, 2)), size = 3.5) +
  facet_wrap(~prefilter) +
  scale_fill_gradient(low = "white", high = eam_palette[1],
                      limits = c(0, 1)) +
  labs(x = "Core Distribution", y = "Contaminant Bounds",
       fill = "F1 Score") +
  clean_plot()

ggsave(here("figures", "10_classification_heatmap.pdf"),
       p_classification, width = 6.5, height = 4, device = "pdf")

## 7.3) Visualization: var_rt recovery (data-driven bounds, with/without prefilter)

# Focus on data-driven bounds (the recommended default) and show the
# effect of pre-filtering across distributions
plot_data_main <- bind_rows(
  sumstats_contaminated |>
    mutate(approach = "none", distribution = "none",
           bound_setting = "none", prefilter = "none"),
  sumstats_all_approaches |>
    filter(bound_setting == "data_driven", approach == "mixture"),
  sumstats_clean |>
    mutate(approach = "ground_truth", distribution = "none",
           bound_setting = "none", prefilter = "none")
) |>
  mutate(
    accuracy = n_upper / n_trials,
    dist_short = str_replace(distribution, "invgaussian", "invGauss") |>
      str_replace("exgaussian", "exGauss") |>
      str_replace("lognormal", "logNorm"),
    label = case_when(
      approach == "none" ~ "Contam.",
      approach == "ground_truth" ~ "Clean",
      prefilter == "none" ~ paste0(dist_short),
      prefilter == "hard" ~ paste0(dist_short, "\n+hard")
    )
  )

# Order labels logically
label_order <- c(
  "Contam.",
  "exGauss", "logNorm", "invGauss",
  "exGauss\n+hard", "logNorm\n+hard", "invGauss\n+hard",
  "Clean"
)
plot_data_main <- plot_data_main |>
  mutate(label = factor(label, levels = label_order))

p_var_rt <- plot_data_main |>
  ggplot(aes(x = label, y = var_rt, fill = label)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  facet_wrap(~condition,
             labeller = labeller(condition = str_to_title)) +
  labs(x = "", y = "Variance of RT") +
  guides(fill = "none") +
  clean_plot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(here("figures", "10_var_rt_by_distribution.pdf"),
       p_var_rt, width = 6.5, height = 5, device = "pdf")

## 7.4) Visualization: all summary stats (data-driven, with/without prefilter) --

p_bias <- plot_data_main |>
  pivot_longer(c(mean_rt, var_rt, accuracy),
               names_to = "statistic", values_to = "value") |>
  mutate(
    statistic = factor(statistic,
                       levels = c("mean_rt", "var_rt", "accuracy"),
                       labels = c("Mean RT", "Var RT", "Accuracy"))
  ) |>
  ggplot(aes(x = label, y = value, fill = label)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  facet_grid(statistic ~ condition, scales = "free_y",
             labeller = labeller(condition = str_to_title)) +
  labs(x = "", y = "") +
  guides(fill = "none") +
  clean_plot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(here("figures", "10_parameter_bias.pdf"),
       p_bias, width = 6.5, height = 7, device = "pdf")

## 7.5) Print recommendation ----------------------------------------------------

cat("\n=== Summary ===\n\n")
cat("Three factors evaluated:\n")
cat("  1. Distribution: exgaussian, lognormal, invgaussian\n")
cat("  2. Contaminant bounds: data-driven, fixed narrow, fixed wide\n")
cat("  3. Pre-filtering: none vs hard (RT < 150ms or RT > 5000ms)\n\n")

cat("CLASSIFICATION (top 5 by F1):\n")
classification_full |>
  select(distribution, bound_setting, prefilter, f1) |>
  arrange(desc(f1)) |>
  mutate(f1 = round(f1, 3)) |>
  slice_head(n = 5) |>
  print()

cat("\nSUMMARY STAT RECOVERY (top 5 by var_rt MAE, mixture method):\n")
bias_summary |>
  filter(approach == "mixture") |>
  select(condition, distribution, bound_setting, prefilter, mae_var_rt) |>
  arrange(condition, mae_var_rt) |>
  group_by(condition) |>
  slice_head(n = 5) |>
  ungroup() |>
  mutate(mae_var_rt = round(mae_var_rt, 4)) |>
  print(n = Inf)
