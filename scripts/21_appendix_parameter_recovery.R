#' Appendix: Parameter Recovery on Simulated Data
#'
#' This script fits the EZ-Diffusion 3par, Censored-Shifted Wald crisk, and DDM to simulated data with
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
source(here("functions", "extract_fixef.R"))

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

# EZ-Diffusion 3par
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
  chains = chains, cores = cores, threads = threading(threads),
  warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_ezdm_3par")
)

# Censored-Shifted Wald crisk
model_cswald <- cswald(rt = "rt", response = "response", version = "crisk")
formula_cswald <- bmf(
  drift ~ 0 + condition + (0 + condition | id),
  bound ~ 1 + (1 | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 1 + (1 | id)
)
fit_cswald_crisk <- bmm(
  formula = formula_cswald, data = data_sim, model = model_cswald,
  chains = chains, cores = cores, threads = threading(threads),
  warmup = warmup, iter = iter,
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
  chains = chains, cores = cores, threads = threading(threads),
  warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_ddm")
)

###############################################################################!
# 1) Group-Level Recovery ------------------------------------------------------
###############################################################################!

# Shared helper (functions/extract_fixef.R): tidy fixef() extraction with
# link-inverse transforms per parameter
group_transforms <- list(bound = exp, ndt = exp)

group_recovery <- bind_rows(
  extract_fixef(fit_ezdm_3par, "EZ-Diffusion (3par)", c("drift", "bound", "ndt"),
                group_transforms),
  extract_fixef(fit_cswald_crisk, "Censored-Shifted Wald (crisk)", c("drift", "bound", "ndt"),
                group_transforms),
  extract_fixef(fit_ddm, "Full DDM", c("drift", "bound", "ndt"), group_transforms)
) |>
  clean_levels() |>
  mutate(level = if_else(level == "", "overall", level))

# Add true values
true_vals <- tibble(
  parameter = c("drift", "drift", "bound", "ndt"),
  level = c("easy", "hard", "overall", "overall"),
  true_value = c(2.5, 1.0, 1.5, 0.3)
)

group_recovery <- group_recovery |>
  left_join(true_vals, by = c("parameter", "level")) |>
  mutate(
    model = factor(model, levels = c("EZ-Diffusion (3par)", "Censored-Shifted Wald (crisk)", "Full DDM")),
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
  # coef()$id carries the 95% CrI bounds (Q2.5/Q97.5) alongside the mean, so
  # we extract them here to report posterior uncertainty per participant.
  drift_easy <- tibble(
    model = model_name,
    parameter = "drift_easy",
    id = 1:n_subjects,
    recovered = co[1:n_subjects, "Estimate", "drift_conditioneasy"],
    lower = co[1:n_subjects, "Q2.5", "drift_conditioneasy"],
    upper = co[1:n_subjects, "Q97.5", "drift_conditioneasy"]
  )

  drift_hard <- tibble(
    model = model_name,
    parameter = "drift_hard",
    id = 1:n_subjects,
    recovered = co[1:n_subjects, "Estimate", "drift_conditionhard"],
    lower = co[1:n_subjects, "Q2.5", "drift_conditionhard"],
    upper = co[1:n_subjects, "Q97.5", "drift_conditionhard"]
  )

  # Bound (log link): exponentiate the point and both CrI bounds. exp() is
  # monotonic, so exp(Q2.5)/exp(Q97.5) are the natural-scale 95% CrI limits.
  bound <- tibble(
    model = model_name,
    parameter = "bound",
    id = 1:n_subjects,
    recovered = exp(co[1:n_subjects, "Estimate", "bound_Intercept"]),
    lower = exp(co[1:n_subjects, "Q2.5", "bound_Intercept"]),
    upper = exp(co[1:n_subjects, "Q97.5", "bound_Intercept"])
  )

  # Ndt (log link): same natural-scale back-transform of the CrI bounds.
  ndt <- tibble(
    model = model_name,
    parameter = "ndt",
    id = 1:n_subjects,
    recovered = exp(co[1:n_subjects, "Estimate", "ndt_Intercept"]),
    lower = exp(co[1:n_subjects, "Q2.5", "ndt_Intercept"]),
    upper = exp(co[1:n_subjects, "Q97.5", "ndt_Intercept"])
  )

  bind_rows(drift_easy, drift_hard, bound, ndt)
}

subject_recovery <- bind_rows(
  extract_subject_params(fit_ezdm_3par, "EZ-Diffusion (3par)"),
  extract_subject_params(fit_cswald_crisk, "Censored-Shifted Wald (crisk)"),
  extract_subject_params(fit_ddm, "Full DDM")
)

# Merge with true subject-level parameters. Drift varies by condition;
# bound_subj and ndt_subj are constant across conditions, so one row per id
# suffices for those.
true_subject <- bind_rows(
  true_params_subject |>
    transmute(id, parameter = paste0("drift_", condition),
              true_value = drift_subj),
  true_params_subject |>
    filter(condition == "easy") |>
    pivot_longer(c(bound_subj, ndt_subj),
                 names_to = "parameter", values_to = "true_value") |>
    transmute(id, parameter = str_remove(parameter, "_subj"), true_value)
)

subject_recovery <- subject_recovery |>
  left_join(true_subject, by = c("id", "parameter")) |>
  mutate(
    model = factor(model, levels = c("EZ-Diffusion (3par)", "Censored-Shifted Wald (crisk)", "Full DDM")),
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
# One facet_wrap row per estimation approach (parameters across columns),
# stacked with patchwork. Each panel keeps free scales so parameters on
# different ranges (notably ndt) stay informative. Model colors follow the
# c(1, 3, 4) palette indices used for the group-level figure.
model_colors <- c("EZ-Diffusion (3par)"    = eam_palette[1],
                  "Censored-Shifted Wald (crisk)" = eam_palette[3],
                  "Full DDM"     = eam_palette[4])

recovery_row <- function(model_name, x_title = FALSE, y_title = FALSE) {
  d    <- filter(subject_recovery, model == model_name)
  cors <- filter(subject_cors, model == model_name)
  ggplot(d, aes(x = true_value, y = recovered)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_pointrange(aes(ymin = lower, ymax = upper),
                    color = model_colors[[model_name]],
                    alpha = 0.5, size = 0.3, linewidth = 0.3) +
    geom_text(data = cors, aes(x = -Inf, y = Inf, label = label),
              color = model_colors[[model_name]],
              hjust = -0.1, vjust = 1.3, size = 3) +
    facet_wrap(~ parameter, nrow = 1, ncol = 4, scales = "free") +
    labs(title = model_name,
         x = if (x_title) "True value" else NULL,
         y = if (y_title) "Recovered estimate" else NULL) +
    clean_plot()
}

p_subject <- recovery_row("EZ-Diffusion (3par)") /
  recovery_row("Censored-Shifted Wald (crisk)", y_title = TRUE) /
  recovery_row("Full DDM", x_title = TRUE)

ggsave(here("figures", "21_recovery_subject.pdf"),
       p_subject, width = 6.5, height = 7.5)
