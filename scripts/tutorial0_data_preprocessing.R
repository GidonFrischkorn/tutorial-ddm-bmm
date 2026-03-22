#' Tutorial 0: The Impact of Contaminants on EAM Parameter Estimates
#'
#' This tutorial demonstrates WHY preprocessing matters for evidence
#' accumulation modeling. Using simulated data with known parameters, we show
#' that contaminant RTs (fast guesses, slow outliers) systematically bias
#' parameter estimates. We then demonstrate bmm's preprocessing tools —
#' flag_contaminant_rts() and ezdm_summary_stats() — and verify that cleaning
#' restores accurate parameter recovery.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# Load both data versions
data_contaminated <- read_csv(here("data", "simulated_rt_data_contaminated.csv"))
data_clean        <- read_csv(here("data", "simulated_rt_data_clean.csv"))

# Load true generating parameters for comparison
true_params <- read_csv(here("data", "true_parameters.csv"))

###############################################################################!
# 1) Bias Demonstration: Contaminants Distort Parameter Estimates --------------
###############################################################################!

# First, compute summary statistics from the contaminated data using the
# simplest method (no contamination handling). This is what a naive analysis
# would produce.
#
# TODO: Implement once ezdm_summary_stats is available
#
# sumstats_contaminated <- ezdm_summary_stats(
#   data_contaminated, rt = rt, response = response,
#   .by = c(id, condition), version = "3par",
#   method = "simple"
# )
#
# # For comparison, compute summary statistics from the clean data
# sumstats_clean <- ezdm_summary_stats(
#   data_clean, rt = rt, response = response,
#   .by = c(id, condition), version = "3par",
#   method = "simple"
# )

## 1.1) Visualize the contaminated RT distributions ----------------------------
#
# The contaminated data should show:
# - A spike of very fast RTs (fast guesses, below typical ndt)
# - A heavy right tail (slow outliers)
# - Lower accuracy in affected conditions (fast guesses at chance)
#
# TODO: RT density plots showing clean vs contaminated distributions
#
# p_rt_density <- data_contaminated |>
#   ggplot(aes(x = rt, fill = factor(is_contaminant))) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~condition) +
#   scale_fill_eam() +
#   labs(x = "Reaction Time (s)", y = "Density",
#        fill = "Contaminant") +
#   clean_plot()

## 1.2) Show how contaminants bias summary statistics --------------------------
#
# Contaminant effects on summary statistics:
# - Fast guesses: decrease mean_rt, INCREASE var_rt dramatically, decrease accuracy
# - Slow outliers: increase mean_rt, INCREASE var_rt
# - Net effect on var_rt is always inflation -> biases EZDM parameters
#
# EZ equations map summary stats to parameters as follows:
# - Inflated var_rt -> overestimated bound (a) and underestimated drift (v)
# - Decreased accuracy -> underestimated drift (v)
# - Biased mean_rt -> biased ndt (Ter)
#
# TODO: Compare summary statistics side-by-side
#
# comparison_sumstats <- bind_rows(
#   sumstats_contaminated |> mutate(source = "Contaminated"),
#   sumstats_clean |> mutate(source = "Clean")
# )
#
# # Plot: mean_rt, var_rt, accuracy by condition and source
# p_sumstats <- comparison_sumstats |>
#   pivot_longer(c(mean_rt, var_rt, n_upper),
#                names_to = "statistic", values_to = "value") |>
#   ggplot(aes(x = source, y = value, fill = source)) +
#   geom_boxplot() +
#   facet_grid(statistic ~ condition, scales = "free_y") +
#   scale_fill_eam() +
#   clean_plot()

###############################################################################!
# 2) Contaminant Detection & Cleaning ------------------------------------------
###############################################################################!

## 2.1) Identifying contaminants with flag_contaminant_rts() -------------------
#
# flag_contaminant_rts() fits a mixture of a core RT distribution (e.g.,
# ex-Gaussian) and a uniform contaminant distribution. It returns a
# contamination probability for each trial.
#
# TODO: Demonstrate flag_contaminant_rts() once the function is merged
#
# data_flagged <- data_contaminated |>
#   group_by(id, condition) |>
#   mutate(
#     contam_prob = flag_contaminant_rts(rt)
#   ) |>
#   ungroup()

## 2.2) Visualize flagged trials -----------------------------------------------
#
# Show RT density with flagged trials highlighted. The function should
# correctly identify the injected contaminants (fast guesses + slow outliers).
#
# TODO: Compare flagged trials against known contaminant labels
#
# # How well does flag_contaminant_rts recover the true contaminants?
# classification <- data_flagged |>
#   mutate(flagged = contam_prob > 0.5) |>
#   group_by(is_contaminant, flagged) |>
#   tally() |>
#   ungroup()
# print(classification)  # confusion matrix

## 2.3) Compute summary statistics with different methods ----------------------
#
# ezdm_summary_stats() offers three methods for handling contamination:
# - "simple": Standard mean/variance (ignores contaminants)
# - "mixture": Fits a mixture model to separate core RT from contaminants
# - "robust": Uses IQR/MAD-based estimators resistant to outliers
#
# The mixture method is recommended because it models the contamination
# process directly, analogous to what flag_contaminant_rts() does.
#
# TODO: Compare methods
#
# data_for_ezdm <- data_flagged |> filter(contam_prob <= 0.5)
#
# sumstats_simple <- ezdm_summary_stats(
#   data_for_ezdm, rt = rt, response = response,
#   .by = c(id, condition), version = "3par", method = "simple"
# )
#
# sumstats_mixture <- ezdm_summary_stats(
#   data_contaminated, rt = rt, response = response,
#   .by = c(id, condition), version = "3par", method = "mixture"
# )
#
# sumstats_robust <- ezdm_summary_stats(
#   data_contaminated, rt = rt, response = response,
#   .by = c(id, condition), version = "3par", method = "robust"
# )

###############################################################################!
# 3) Bias Reduction Verification -----------------------------------------------
###############################################################################!

# Compare parameter estimates from contaminated vs cleaned data against
# the known generating values. This is the payoff: preprocessing visibly
# reduces bias and brings estimates closer to ground truth.
#
# TODO: Fit EZDM to contaminated and cleaned summary stats, compare
#
# # Quick EZ-equation parameter estimates (no Bayesian fitting needed here)
# # The EZ equations provide closed-form parameter estimates from summary stats.
# # This allows a fast comparison without full model fitting.
#
# compare_params <- bind_rows(
#   # From contaminated data (naive)
#   sumstats_contaminated |> mutate(source = "Contaminated (naive)"),
#   # From cleaned data (after flag_contaminant_rts + simple aggregation)
#   sumstats_simple |> mutate(source = "Cleaned (flag + simple)"),
#   # From contaminated data with mixture method (no pre-filtering needed)
#   sumstats_mixture |> mutate(source = "Mixture method (no filtering)")
# )

## 3.1) Parameter bias plot ----------------------------------------------------
#
# Figure: Three panels (drift, bound, ndt), each showing the distribution of
# individual-level parameter estimates across participants for each method,
# with the true generating value as a reference line.
#
# This figure is the key takeaway: contaminants bias all parameters, and
# proper preprocessing (via flag_contaminant_rts or the mixture method in
# ezdm_summary_stats) substantially reduces this bias.
#
# TODO: Create bias comparison figure
#
# p_bias <- compare_params |>
#   pivot_longer(c(drift_est, bound_est, ndt_est),
#                names_to = "parameter", values_to = "estimate") |>
#   left_join(true_params_long, by = c("condition", "parameter")) |>
#   ggplot(aes(x = source, y = estimate - true_value)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#   geom_boxplot(aes(fill = source)) +
#   facet_grid(parameter ~ condition, scales = "free_y") +
#   scale_fill_eam() +
#   labs(x = "", y = "Bias (estimate - true)") +
#   clean_plot(axis.text.x = element_text(angle = 30, hjust = 1))

## 3.2) Summary ----------------------------------------------------------------
#
# Key takeaways for the reader:
# 1. Contaminant RTs (even ~7%) produce substantial parameter bias
# 2. Fast guesses primarily inflate var_rt -> bias drift and bound
# 3. flag_contaminant_rts() effectively identifies contaminant trials
# 4. The mixture method in ezdm_summary_stats() handles contamination
#    without requiring a separate filtering step
# 5. After cleaning, parameter estimates closely match ground truth
