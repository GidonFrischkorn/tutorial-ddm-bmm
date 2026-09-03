#' Tutorial 0: The Impact of Contaminants on EAM Parameter Estimates
#'
#' This tutorial demonstrates WHY preprocessing matters for evidence
#' accumulation modeling. Using simulated data with known parameters, we show
#' that contaminant RTs systematically bias the summary statistics that EAMs
#' rely on, and that bmm's mixture-based preprocessing tools --
#' flag_contaminant_rts() and ezdm_summary_stats() -- restore accurate summary
#' statistics. The contaminants come from the three processes evaluated in the
#' Monte Carlo study (scripts/02_simulation_contaminant_methods.R):
#' leading-edge anticipations, delayed start-ups, and informationless
#' responses.
#'
#' The systematic comparison of flagging settings (core distributions,
#' contaminant bounds, hard pre-filtering, threshold vs. probabilistic
#' exclusion) lives in scripts/10b_preprocessing_comparison.R.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, tidyverse, patchwork)
source(here("functions", "clean_plot.R"))

# Simulated data from scripts/20_prepare_simulated_data.R:
# 50 participants x 2 conditions x 200 trials with known DDM parameters.
# The contaminated version injects 3% of trials from each of three processes
# (9% in total); is_contaminant flags the injected trials and
# contaminant_type records which process produced them (ground truth).
data_contaminated <- read_csv(
  here("data", "simulated_rt_data_contaminated.csv"),
  show_col_types = FALSE
)
data_clean <- read_csv(
  here("data", "simulated_rt_data_clean.csv"),
  show_col_types = FALSE
)

###############################################################################!
# 1) Bias Demonstration: Contaminants Distort Summary Statistics ---------------
###############################################################################!

## 1.1) Visualize where each contaminant process falls -------------------------

# The three processes differ in how far they fall from the RT distribution the
# decision process produces, which is what determines whether a preprocessing
# method can separate them (Appendix A). Each panel shows the shape of the core
# distribution and the shape of the injected contaminants; both are scaled to
# the same area, so the panels compare where the components sit, not how many
# trials each contains (contaminants are 3% of trials).

process_levels <- c("leading-edge anticipation", "delayed start-up",
                    "informationless")
process_labels <- c("Leading-edge anticipations", "Delayed start-ups",
                    "Informationless responses")

components <- data_contaminated |>
  filter(contaminant_type != "none") |>
  transmute(condition, rt,
            component = factor(contaminant_type, levels = process_levels,
                               labels = process_labels)) |>
  bind_rows(
    # the core distribution is the same in every panel
    expand_grid(data_clean |> select(condition, rt),
                component = factor(process_labels, levels = process_labels)) |>
      mutate(source = "Decision process")
  ) |>
  mutate(
    source = if_else(is.na(source), "Contaminants", source),
    source = factor(source, levels = c("Decision process", "Contaminants")),
    condition = factor(condition, levels = c("easy", "hard"),
                       labels = c("Easy", "Hard"))
  )

p_rt_density <- components |>
  ggplot(aes(x = rt, fill = source, color = source)) +
  geom_density(alpha = 0.55, linewidth = 0.4) +
  facet_grid(condition ~ component, scales = "free_y") +
  scale_fill_manual(values = c("Decision process" = "grey65",
                               "Contaminants" = eam_palette[4])) +
  scale_color_manual(values = c("Decision process" = "grey35",
                                "Contaminants" = eam_palette[4])) +
  coord_cartesian(xlim = c(0, 3.5)) +
  labs(x = "Reaction Time (s)", y = "Density", fill = NULL, color = NULL) +
  clean_plot() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(here("figures", "10_rt_density.pdf"),
       p_rt_density, width = 6.5, height = 4, device = "pdf")

## 1.2) Compare summary statistics: contaminated vs clean ----------------------

# Compute EZDM summary stats from contaminated data (no contaminant handling)
sumstats_contaminated <- data_contaminated |>
  summarize(
    ezdm_summary_stats(rt, response,
                       version = "3par", method = "simple"),
    .by = c(id, condition)
  ) |>
  mutate(source = "Contaminated")

# Compute EZDM summary stats from clean data (ground truth)
sumstats_clean <- data_clean |>
  summarize(
    ezdm_summary_stats(rt, response,
                       version = "3par", method = "simple"),
    .by = c(id, condition)
  ) |>
  mutate(source = "Clean")

# Show how contaminants distort the summary statistics that EAMs
# rely on: mean_rt, var_rt, and accuracy (n_upper / n_trials)
sumstats_comparison <- bind_rows(
  sumstats_contaminated, sumstats_clean
) |>
  mutate(accuracy = n_upper / n_trials)

p_sumstats <- sumstats_comparison |>
  pivot_longer(c(mean_rt, var_rt, accuracy),
               names_to = "statistic", values_to = "value") |>
  mutate(
    statistic = factor(statistic,
                       levels = c("mean_rt", "var_rt", "accuracy"),
                       labels = c("Mean RT", "Var RT", "Accuracy"))
  ) |>
  ggplot(aes(x = source, y = value, fill = source)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  facet_grid(statistic ~ condition, scales = "free_y",
             labeller = labeller(condition = str_to_title)) +
  scale_fill_manual(values = eam_palette[c(4, 1)]) +
  labs(x = "", y = "") +
  guides(fill = "none") +
  clean_plot()

ggsave(here("figures", "10_sumstats_bias.pdf"),
       p_sumstats, width = 6.5, height = 5, device = "pdf")

###############################################################################!
# 2) Contaminant Detection with Default Settings --------------------------------
###############################################################################!

# flag_contaminant_rts() decomposes each participant x condition RT
# distribution into a core component (here: ex-Gaussian) and a uniform
# contaminant component, returning the probability that each trial is a
# contaminant. contaminant_bound = c("min", "max") expects contaminants at
# both ends of the distribution (premature responses AND slow outliers).
data_flagged <- data_contaminated |>
  group_by(id, condition) |>
  mutate(
    contam_prob = flag_contaminant_rts(
      rt,
      distribution = "exgaussian",
      contaminant_bound = c("min", "max"),
      maxit = 500
    )
  ) |>
  ungroup() |>
  mutate(flagged = contam_prob > 0.5)

# Because the data are simulated, we can score the classification against
# the known ground truth (is_contaminant)
classification <- data_flagged |>
  summarize(
    sensitivity = sum(is_contaminant & flagged, na.rm = TRUE) /
      sum(is_contaminant),
    specificity = sum(!is_contaminant & !flagged, na.rm = TRUE) /
      sum(!is_contaminant),
    precision   = sum(is_contaminant & flagged, na.rm = TRUE) /
      sum(flagged, na.rm = TRUE)
  )
cat("\n=== Default flagging vs. ground truth ===\n")
cat("Sensitivity (contaminants caught):     ",
    round(classification$sensitivity, 3), "\n")
cat("Specificity (clean trials kept):       ",
    round(classification$specificity, 3), "\n")
cat("Precision (flagged = true contaminant):",
    round(classification$precision, 3), "\n")

# The scores show the asymmetry the manuscript discusses: flags are rare among
# clean trials (high specificity), but contaminants overlapping the core --
# leading-edge anticipations and informationless responses -- largely escape
# detection. validate_fast_guesses() offers an additional check that flagged
# fast trials show chance accuracy; see scripts/10b_preprocessing_comparison.R
# for it and for the systematic comparison of flagging settings.

###############################################################################!
# 3) Cleaning Restores Accurate Summary Statistics -----------------------------
###############################################################################!

# The mixture method in ezdm_summary_stats() integrates contaminant handling
# and aggregation in one step: probable contaminants are down-weighted rather
# than excluded with a hard cutoff.
sumstats_mixture <- data_contaminated |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par",
                       method = "mixture", distribution = "exgaussian",
                       maxit = 500),
    .by = c(id, condition)
  ) |>
  mutate(source = "Mixture-cleaned")

# Median absolute error relative to the clean ground truth
mae_summary <- bind_rows(
  sumstats_contaminated |> mutate(source = "Contaminated (no handling)"),
  sumstats_mixture
) |>
  left_join(
    sumstats_clean |>
      select(id, condition, mean_rt_clean = mean_rt, var_rt_clean = var_rt),
    by = c("id", "condition")
  ) |>
  summarize(
    mae_mean_rt = median(abs(mean_rt - mean_rt_clean)),
    mae_var_rt  = median(abs(var_rt - var_rt_clean)),
    .by = c(source, condition)
  )

cat("\n=== Summary statistic recovery (median absolute error) ===\n")
mae_summary |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(condition, mae_var_rt) |>
  print()

# The mixture method substantially reduces the error that contaminants
# introduce in mean RT and RT variance -- and by extension the bias in the
# EZ-Diffusion parameters computed from them.
