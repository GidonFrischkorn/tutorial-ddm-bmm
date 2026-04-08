#' Prepare Wagenmakers et al. (2008) Lexical Decision Data
#'
#' Shared data preparation script sourced by all tutorial scripts.
#' Loads the speed_acc dataset from rtdists, creates accuracy-coded and
#' response-coded versions, applies contaminant detection, and computes
#' EZDM summary statistics.
#'
#' After sourcing this script, the following objects are available:
#'   - data_lexical:        Raw data (factors converted, errors/censored removed)
#'   - data_accuracy:       Accuracy-coded (upper = correct, lower = error)
#'   - data_response:       Response-coded (upper = "word", lower = "nonword")
#'   - data_accuracy_clean: Accuracy-coded, contaminants removed
#'   - data_response_clean: Response-coded, contaminants removed
#'   - data_acc_ezdm:       EZDM 3par summary stats (accuracy-coded)
#'   - data_resp_ezdm:      EZDM 4par summary stats (response-coded)

# -- Load raw data ------------------------------------------------------------
data("speed_acc", package = "rtdists")

data_lexical <- speed_acc |>
  mutate(
    response = as.character(response),
    stim_cat = as.character(stim_cat),
    condition = factor(condition, levels = c("speed", "accuracy")),
    # Create 4-level frequency variable matching Wagenmakers et al. (2008)
    # The raw variable has 6 levels (3 word + 3 nonword-derived-from).
    # Collapse all nonword frequencies into a single "nonword" category.
    freq = factor(
      ifelse(stim_cat == "nonword", "nonword", as.character(frequency)),
      levels = c("high", "low", "very_low", "nonword")
    )
  ) |>
  filter(response != "error", censor == FALSE)

# -- Validate mixture model against published hard cutoffs --------------------
# The original study excluded RT < 180ms and RT > 3000ms (the censor variable).
# We apply flag_contaminant_rts() to ALL data (including censored trials)
# to compare the two approaches.

data_all <- speed_acc |>
  mutate(
    response = as.character(response),
    stim_cat = as.character(stim_cat),
    condition = factor(condition, levels = c("speed", "accuracy")),
    freq = factor(
      ifelse(stim_cat == "nonword", "nonword", as.character(frequency)),
      levels = c("high", "low", "very_low", "nonword")
    )
  ) |>
  filter(response != "error")

data_validation <- data_all |>
  group_by(id, condition, stim_cat) |>
  mutate(
    contam_prob = flag_contaminant_rts(
      rt, distribution = "exgaussian",
      contaminant_bound = c("min", "max")
    )
  ) |>
  ungroup() |>
  mutate(
    flagged_mixture = contam_prob > 0.5,
    flagged_censor = censor
  )

# Confusion matrix: censor (hard bounds) vs mixture model
cat("\n=== Censor vs Mixture Model Comparison ===\n")
cat("Both flag:    ",
    sum(data_validation$flagged_censor & data_validation$flagged_mixture), "\n")
cat("Censor only:  ",
    sum(data_validation$flagged_censor & !data_validation$flagged_mixture), "\n")
cat("Mixture only: ",
    sum(!data_validation$flagged_censor & data_validation$flagged_mixture), "\n")
cat("Neither:      ",
    sum(!data_validation$flagged_censor & !data_validation$flagged_mixture), "\n")

tp <- sum(data_validation$flagged_censor & data_validation$flagged_mixture)
fn <- sum(data_validation$flagged_censor & !data_validation$flagged_mixture)
cat("Sensitivity:  ", round(tp / (tp + fn), 3), "\n\n")

rm(data_all, data_validation)

# -- Filter to non-censored, non-error trials for the tutorials ---------------
# (From here on, we use the standard filtered data)

# -- Two response codings from the same data ----------------------------------

# Accuracy-coded: correct = upper, error = lower
data_accuracy <- data_lexical |>
  mutate(response = ifelse(response == stim_cat, "upper", "lower"))

# Response-coded: "word" response = upper, "nonword" response = lower
data_response <- data_lexical |>
  mutate(response = ifelse(response == "word", "upper", "lower"))

# -- Contaminant detection ----------------------------------------------------
data_accuracy <- data_accuracy |>
  group_by(id, condition, stim_cat) |>
  mutate(
    contam_prob = flag_contaminant_rts(
      rt,
      distribution = "exgaussian",
      contaminant_bound = c("min", "max")
    )
  ) |>
  ungroup()

# Clean trial-level data (contaminants removed)
# Two approaches are available for trial exclusion:
#   Threshold-based: exclude if contam_prob > 0.5 (deterministic)
#   Probabilistic:   exclude with probability = contam_prob (stochastic)
# We use threshold-based as the default for reproducibility.
# See 10_tutorial0_data_preprocessing.R for a comparison of both approaches.
data_accuracy_clean <- data_accuracy |> filter(contam_prob <= 0.5)
data_response_clean <- data_response |>
  semi_join(data_accuracy_clean, by = c("id", "condition", "stim_cat", "rt"))

# -- EZDM summary statistics --------------------------------------------------

# 3par (accuracy-coded)
data_acc_ezdm <- data_accuracy |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par",
                       method = "mixture", distribution = "exgaussian"),
    .by = c(id, condition, stim_cat)
  )

# 4par (response-coded)
data_resp_ezdm <- data_response |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par",
                       method = "mixture", distribution = "exgaussian"),
    .by = c(id, condition, stim_cat)
  )

# -- Save prepared datasets to data/ ------------------------------------------
# Trial-level datasets include freq for use in DDM tutorials
write_csv(data_accuracy_clean, here("data", "lexical_accuracy_triallevel.csv"))
write_csv(data_response_clean, here("data", "lexical_response_triallevel.csv"))
write_csv(data_acc_ezdm, here("data", "lexical_accuracy_ezdm3par.csv"))
write_csv(data_resp_ezdm, here("data", "lexical_response_ezdm4par.csv"))
