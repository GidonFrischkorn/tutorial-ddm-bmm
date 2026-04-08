#' Appendix: Parameter Recovery on Simulated Data
#'
#' This script fits EZDM 3par, CSWALD crisk, and DDM to simulated data with
#' known generating parameters, then evaluates recovery at both the group level
#' (population means) and subject level (individual differences).
#'
#' Requires: Run scripts/20_prepare_simulated_data.R first to generate the
#' simulated datasets in data/.

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
warmup <- 1000
iter   <- 2000
adapt_delta   <- 0.90
max_treedepth <- 12

# -- Load simulated data ------------------------------------------------------
data_sim <- read_csv(here("data", "simulated_rt_data_clean.csv"),
                     show_col_types = FALSE) |>
  mutate(response = ifelse(response == 1, "upper", "lower"))

data_ezdm <- read_csv(here("data", "simulated_ezdm_3par.csv"),
                       show_col_types = FALSE)

# -- Load true generating parameters -----------------------------------------
true_params <- read_csv(here("data", "true_parameters.csv"),
                         show_col_types = FALSE)
true_params_subject <- read_csv(here("data", "true_parameters_subject.csv"),
                                 show_col_types = FALSE)

###############################################################################!
# 0b) Fit Models ---------------------------------------------------------------
###############################################################################!

# EZDM 3par
model_ezdm <- ezdm(mean_rt = "mean_rt", var_rt = "var_rt",
                    n_upper = "n_upper", n_trials = "n_trials",
                    version = "3par")
formula_ezdm <- bmf(
  drift ~ 0 + condition + (0 + condition | id),
  bound ~ 1 + (1 | id),
  ndt   ~ 1 + (1 | id)
)
fit_ezdm_3par <- bmm(
  formula = formula_ezdm, data = data_ezdm, model = model_ezdm,
  chains = chains, cores = cores, warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_ezdm_3par")
)

# CSWALD crisk
model_cswald <- cswald(rt = "rt", response = "response", version = "crisk")
formula_cswald <- bmf(
  drift ~ 0 + condition + (0 + condition | id),
  bound ~ 1 + (1 | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 1 + (1 | id)
)
fit_cswald_crisk <- bmm(
  formula = formula_cswald, data = data_sim, model = model_cswald,
  chains = chains, cores = cores, warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_cswald_crisk")
)

# DDM
model_ddm <- ddm(rt = "rt", response = "response")
formula_ddm <- bmf(
  drift ~ 0 + condition + (0 + condition | id),
  bound ~ 1 + (1 | id),
  ndt   ~ 1 + (1 | id)
)
fit_ddm <- bmm(
  formula = formula_ddm, data = data_sim, model = model_ddm,
  chains = chains, cores = cores, warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_ddm")
)

###############################################################################!
# 1) Group-Level Recovery ------------------------------------------------------
###############################################################################!

extract_group_params <- function(fit, model_name) {
  fe <- fixef(fit)

  # Drift (identity link)
  drift_rows <- grep("^drift_", rownames(fe), value = TRUE)
  drift_est <- tibble(
    model = model_name,
    parameter = "drift",
    level = sub("^drift_condition", "", drift_rows),
    estimate = fe[drift_rows, "Estimate"],
    lower = fe[drift_rows, "Q2.5"],
    upper = fe[drift_rows, "Q97.5"]
  )

  # Bound (log link)
  bound_row <- grep("^bound_", rownames(fe), value = TRUE)
  bound_est <- tibble(
    model = model_name,
    parameter = "bound",
    level = "overall",
    estimate = exp(fe[bound_row, "Estimate"]),
    lower = exp(fe[bound_row, "Q2.5"]),
    upper = exp(fe[bound_row, "Q97.5"])
  )

  # Ndt (log link)
  ndt_row <- grep("^ndt_", rownames(fe), value = TRUE)
  ndt_est <- tibble(
    model = model_name,
    parameter = "ndt",
    level = "overall",
    estimate = exp(fe[ndt_row, "Estimate"]),
    lower = exp(fe[ndt_row, "Q2.5"]),
    upper = exp(fe[ndt_row, "Q97.5"])
  )

  bind_rows(drift_est, bound_est, ndt_est)
}

group_recovery <- bind_rows(
  extract_group_params(fit_ezdm_3par, "EZDM 3par"),
  extract_group_params(fit_cswald_crisk, "CSWALD crisk"),
  extract_group_params(fit_ddm, "DDM")
)

# Add true values
true_vals <- tibble(
  parameter = c("drift", "drift", "bound", "ndt"),
  level = c("easy", "hard", "overall", "overall"),
  true_value = c(2.5, 1.0, 1.5, 0.3)
)

group_recovery <- group_recovery |>
  left_join(true_vals, by = c("parameter", "level")) |>
  mutate(
    model = factor(model, levels = c("EZDM 3par", "CSWALD crisk", "DDM")),
    parameter = factor(parameter, levels = c("drift", "bound", "ndt"))
  )

# -- Plot group-level recovery ------------------------------------------------
p_group <- ggplot(group_recovery, aes(x = level, y = estimate,
                                       color = model, shape = model)) +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", alpha = 0.5) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  facet_wrap(~ parameter, scales = "free",
             labeller = labeller(
               parameter = c(drift = "Drift",
                              bound = "Bound",
                              ndt = "NDT"))) +
  scale_color_manual(values = eam_palette[c(1, 3, 4)]) +
  labs(x = NULL, y = "Estimate", color = "Model", shape = "Model") +
  clean_plot()

ggsave(here("figures", "21_recovery_group.pdf"),
       p_group, width = 6.5, height = 3.5)

###############################################################################!
# 2) Subject-Level Recovery ----------------------------------------------------
###############################################################################!

extract_subject_params <- function(fit, model_name, n_subjects = 50) {
  co <- coef(fit)$id

  # With (0 + condition | id) on drift, coef() returns condition-specific
  # subject-level estimates directly.
  drift_easy <- tibble(
    model = model_name,
    parameter = "drift_easy",
    id = 1:n_subjects,
    recovered = co[1:n_subjects, "Estimate", "drift_conditioneasy"]
  )

  drift_hard <- tibble(
    model = model_name,
    parameter = "drift_hard",
    id = 1:n_subjects,
    recovered = co[1:n_subjects, "Estimate", "drift_conditionhard"]
  )

  # Bound (log link)
  bound <- tibble(
    model = model_name,
    parameter = "bound",
    id = 1:n_subjects,
    recovered = exp(co[1:n_subjects, "Estimate", "bound_Intercept"])
  )

  # Ndt (log link)
  ndt <- tibble(
    model = model_name,
    parameter = "ndt",
    id = 1:n_subjects,
    recovered = exp(co[1:n_subjects, "Estimate", "ndt_Intercept"])
  )

  bind_rows(drift_easy, drift_hard, bound, ndt)
}

subject_recovery <- bind_rows(
  extract_subject_params(fit_ezdm_3par, "EZDM 3par"),
  extract_subject_params(fit_cswald_crisk, "CSWALD crisk"),
  extract_subject_params(fit_ddm, "DDM")
)

# Merge with true subject-level parameters
# Data is long format: id, condition, drift_subj, bound_subj, ndt_subj
true_subject <- bind_rows(
  true_params_subject |>
    filter(condition == "easy") |>
    transmute(id, parameter = "drift_easy", true_value = drift_subj),
  true_params_subject |>
    filter(condition == "hard") |>
    transmute(id, parameter = "drift_hard", true_value = drift_subj),
  true_params_subject |>
    filter(condition == "easy") |>
    transmute(id, parameter = "bound", true_value = bound_subj),
  true_params_subject |>
    filter(condition == "easy") |>
    transmute(id, parameter = "ndt", true_value = ndt_subj)
)

subject_recovery <- subject_recovery |>
  left_join(true_subject, by = c("id", "parameter")) |>
  mutate(
    model = factor(model, levels = c("EZDM 3par", "CSWALD crisk", "DDM")),
    parameter = factor(parameter,
                       levels = c("drift_easy", "drift_hard", "bound", "ndt"),
                       labels = c("Drift (easy)", "Drift (hard)", "Bound", "NDT"))
  )

# Compute correlations
subject_cors <- subject_recovery |>
  group_by(model, parameter) |>
  summarise(r = cor(true_value, recovered), .groups = "drop") |>
  mutate(label = sprintf("r = %.2f", r))

# -- Plot subject-level recovery -----------------------------------------------
p_subject <- ggplot(subject_recovery,
                     aes(x = true_value, y = recovered, color = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text(data = subject_cors,
            aes(x = -Inf, y = Inf, label = label, color = model),
            hjust = -0.1, vjust = 1.3, size = 3, show.legend = FALSE) +
  facet_grid(model ~ parameter, scales = "free") +
  scale_color_manual(values = eam_palette[c(1, 3, 4)]) +
  labs(x = "True value", y = "Recovered estimate") +
  clean_plot() +
  theme(legend.position = "none")

ggsave(here("figures", "21_recovery_subject.pdf"),
       p_subject, width = 6.5, height = 6)
