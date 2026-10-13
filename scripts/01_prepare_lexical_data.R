# Data preparation: Wagenmakers et al. (2008) lexical decision data
# Loads speed_acc (rtdists), creates accuracy- and response-coded versions,
# removes contaminants, and computes the EZ-Diffusion summary statistics.
# Output: trial-level and summary-statistic CSVs in data/.

pacman::p_load(here, bmm, tidyverse, rtdists)

# -- Load raw data ------------------------------------------------------------
data("speed_acc", package = "rtdists")

# Recording failures (response == "error") and the trials the original study
# censored are removed; nonwords are collapsed into one frequency level.
data_lexical <- speed_acc |>
  mutate(
    response = as.character(response),
    stim_cat = as.character(stim_cat),
    condition = factor(condition, levels = c("speed", "accuracy")),
    freq = factor(
      ifelse(stim_cat == "nonword", "nonword", as.character(frequency)),
      levels = c("high", "low", "very_low", "nonword")
    )
  ) |>
  filter(response != "error", censor == FALSE)

# -- Two response codings of the same data -------------------------------------

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

# NA marks a cell whose EM did not converge; the filter would drop it silently.
# Probabilistic alternative (set a seed first; needs a unique trial id for the join):
#   data_accuracy_clean <- data_accuracy |> filter(runif(n()) > contam_prob)
stopifnot(!anyNA(data_accuracy$contam_prob))
data_accuracy_clean <- data_accuracy |> filter(contam_prob <= 0.5)
data_response_clean <- data_response |>
  semi_join(data_accuracy_clean, by = c("id", "condition", "stim_cat", "rt"))

# -- EZ-Diffusion summary statistics ------------------------------------------

# 3par (accuracy-coded)
data_acc_ezdm <- data_accuracy_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par",
                       method = "mixture", distribution = "exgaussian"),
    .by = c(id, condition, stim_cat)
  )

# 4par (response-coded): robust method, because some cells have too few
# minority-boundary responses to fit a mixture
data_resp_ezdm <- data_response |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par",
                       method = "robust"),
    .by = c(id, condition, stim_cat)
  )

# -- Save prepared datasets ---------------------------------------------------
write_csv(data_accuracy_clean, here("data", "lexical_accuracy_triallevel.csv"))
write_csv(data_response_clean, here("data", "lexical_response_triallevel.csv"))
write_csv(data_acc_ezdm, here("data", "lexical_accuracy_ezdm3par.csv"))
write_csv(data_resp_ezdm, here("data", "lexical_response_ezdm4par.csv"))
