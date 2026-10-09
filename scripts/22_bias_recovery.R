#' Part I: Parameter Recovery With a Biased Starting Point
#'
#' Fits the 4-parameter EZ-Diffusion Model, the 3-parameter EZ-Diffusion Model
#' (starting point fixed at 0.5), and the full DDM to simulated data whose
#' relative starting point is 0.60 or 0.40 by bias condition, then evaluates
#' recovery at the group level (population means) and the subject level.
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
data_bias <- read_csv(here("data", "simulated_rt_data_bias.csv"),
                      show_col_types = FALSE) |>
  mutate(stimulus = factor(stimulus, levels = c("A", "B")),
         bias     = factor(bias, levels = c("A", "B")))

true_bias <- read_csv(here("data", "true_parameters_bias.csv"),
                      show_col_types = FALSE)
true_bias_subject <- read_csv(here("data", "true_parameters_bias_subject.csv"),
                              show_col_types = FALSE)

# -- Summary statistics for the EZ-Diffusion fits ------------------------------
# One cell per participant x stimulus x bias condition. The 4par version needs
# boundary-specific moments; the 3par version pools both boundaries.
data_bias_4par <- data_bias |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par",
                       method = "mixture", distribution = "exgaussian",
                       maxit = 500),
    .by = c(id, stimulus, bias)
  )
data_bias_3par <- data_bias |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par",
                       method = "mixture", distribution = "exgaussian",
                       maxit = 500),
    .by = c(id, stimulus, bias)
  )
stopifnot(!anyNA(data_bias_4par$contaminant_prop_upper),
          !anyNA(data_bias_4par$contaminant_prop_lower),
          !anyNA(data_bias_3par$contaminant_prop))

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

model_ezdm_4par <- ezdm(
  mean_rt  = c("mean_rt_upper", "mean_rt_lower"),
  var_rt   = c("var_rt_upper", "var_rt_lower"),
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "4par"
)
model_ezdm_3par <- ezdm(
  mean_rt  = "mean_rt",
  var_rt   = "var_rt",
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "3par"
)
model_ddm <- ddm(rt = "rt", response = "response")

# Stimulus type on drift, bias condition on zr; bound and ndt do not vary by
# condition in the generating model.
formula_bias <- bmf(
  drift ~ 0 + stimulus + (0 + stimulus | id),
  bound ~ 1 + (1 | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + bias + (0 + bias | id)
)
# The 3par version has no zr: the same formula without that line.
formula_bias_3par <- bmf(
  drift ~ 0 + stimulus + (0 + stimulus | id),
  bound ~ 1 + (1 | id),
  ndt   ~ 1 + (1 | id)
)

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

fit_ezdm_4par_bias <- bmm(
  formula = formula_bias, data = data_bias_4par, model = model_ezdm_4par,
  chains = chains, cores = cores, threads = threading(threads),
  warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_ezdm_4par_bias")
)

fit_ezdm_3par_bias <- bmm(
  formula = formula_bias_3par, data = data_bias_3par, model = model_ezdm_3par,
  chains = chains, cores = cores, threads = threading(threads),
  warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_ezdm_3par_bias")
)

fit_ddm_bias <- bmm(
  formula = formula_bias, data = data_bias, model = model_ddm,
  chains = chains, cores = cores, threads = threading(threads),
  warmup = warmup, iter = iter,
  backend = "cmdstanr", sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file = here("output", "fit_ddm_bias")
)

fits <- list(
  "EZ-Diffusion (4par)" = fit_ezdm_4par_bias,
  "EZ-Diffusion (3par)" = fit_ezdm_3par_bias,
  "Full DDM"            = fit_ddm_bias
)
model_levels <- names(fits)

# Convergence and wall-clock (slowest chain), the statistics the text reports
for (m in model_levels) {
  cat(sprintf("%-22s max Rhat = %.3f   wall-clock = %.1f min\n", m,
              max(brms::rhat(fits[[m]]), na.rm = TRUE),
              max(rowSums(rstan::get_elapsed_time(fits[[m]]$fit))) / 60))
}

###############################################################################!
# 3) Group-Level Recovery ------------------------------------------------------
###############################################################################!

group_transforms <- list(bound = exp, ndt = exp, zr = plogis)

group_recovery <- bind_rows(
  extract_fixef(fit_ezdm_4par_bias, "EZ-Diffusion (4par)",
                c("drift", "bound", "ndt", "zr"), group_transforms),
  extract_fixef(fit_ezdm_3par_bias, "EZ-Diffusion (3par)",
                c("drift", "bound", "ndt"), group_transforms),
  extract_fixef(fit_ddm_bias, "Full DDM",
                c("drift", "bound", "ndt", "zr"), group_transforms)
) |>
  clean_levels() |>
  # clean_levels() knows the lexical factors; strip this design's own prefixes
  mutate(level = str_remove(level, "^(stimulus|bias)"),
         level = if_else(level == "", "overall", level))

# Generating values: drift by stimulus, zr by bias condition, bound and ndt
# constant. The 3par model has no zr row and is drawn without one.
true_vals <- bind_rows(
  true_bias |> distinct(level = stimulus, true_value = drift) |>
    mutate(parameter = "drift"),
  true_bias |> distinct(level = bias, true_value = zr) |>
    mutate(parameter = "zr"),
  tibble(parameter = c("bound", "ndt"), level = "overall",
         true_value = c(unique(true_bias$bound), unique(true_bias$ndt)))
)

group_recovery <- group_recovery |>
  left_join(true_vals, by = c("parameter", "level")) |>
  mutate(
    model = factor(model, levels = model_levels),
    parameter = factor(parameter, levels = c("drift", "bound", "ndt", "zr"))
  )
stopifnot(!anyNA(group_recovery$true_value))

p_group <- ggplot(group_recovery, aes(x = level, y = estimate,
                                      color = model, shape = model)) +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", alpha = 0.5) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  facet_wrap(~ parameter, scales = "free", nrow = 1,
             labeller = labeller(
               parameter = c(drift = "Drift (stimulus)",
                             bound = "Bound",
                             ndt = "NDT",
                             zr = "Start (bias)"))) +
  scale_color_manual(values = eam_palette[c(1, 2, 4)]) +
  labs(x = NULL, y = "Estimate", color = "Model", shape = "Model") +
  clean_plot()

ggsave(here("figures", "22_bias_recovery_group.pdf"),
       p_group, width = 6.5, height = 3.2)

###############################################################################!
# 4) Subject-Level Recovery ----------------------------------------------------
###############################################################################!

# coef()$id returns condition-specific subject-level estimates with their
# 95% CrI bounds; back-transform log-link (bound, ndt) and logit-link (zr)
# parameters, point and interval alike.
extract_subject_params <- function(fit, model_name, n_subjects = 50) {
  co <- coef(fit)$id
  pars <- dimnames(co)[[3]]
  one <- function(par, label, f = identity) {
    tibble(model = model_name, parameter = label, id = 1:n_subjects,
           recovered = f(co[1:n_subjects, "Estimate", par]),
           lower = f(co[1:n_subjects, "Q2.5", par]),
           upper = f(co[1:n_subjects, "Q97.5", par]))
  }
  out <- bind_rows(
    one("drift_stimulusA", "drift_A"),
    one("drift_stimulusB", "drift_B"),
    one("bound_Intercept", "bound", exp),
    one("ndt_Intercept",   "ndt",   exp)
  )
  if ("zr_biasA" %in% pars) {
    out <- bind_rows(out,
      one("zr_biasA", "zr_A", plogis),
      one("zr_biasB", "zr_B", plogis))
  }
  out
}

subject_recovery <- purrr::map(model_levels, \(m) extract_subject_params(fits[[m]], m)) |>
  bind_rows()

true_subject <- bind_rows(
  true_bias_subject |>
    filter(bias == "A") |>
    transmute(id, parameter = paste0("drift_", stimulus), true_value = drift_subj),
  true_bias_subject |>
    filter(stimulus == "A") |>
    transmute(id, parameter = paste0("zr_", bias), true_value = zr_subj),
  true_bias_subject |>
    filter(stimulus == "A", bias == "A") |>
    pivot_longer(c(bound_subj, ndt_subj),
                 names_to = "parameter", values_to = "true_value") |>
    transmute(id, parameter = str_remove(parameter, "_subj"), true_value)
)

param_levels <- c("drift_A", "drift_B", "bound", "ndt", "zr_A", "zr_B")
# short enough for six panels across 6.5 in (the caption spells them out)
param_labels <- c("Drift (A)", "Drift (B)", "Bound", "NDT",
                  "zr (A)", "zr (B)")

subject_recovery <- subject_recovery |>
  left_join(true_subject, by = c("id", "parameter")) |>
  mutate(model = factor(model, levels = model_levels),
         parameter = factor(parameter, levels = param_levels, labels = param_labels))
stopifnot(!anyNA(subject_recovery$true_value))

subject_cors <- subject_recovery |>
  group_by(model, parameter) |>
  summarise(r = cor(true_value, recovered), .groups = "drop") |>
  mutate(label = sprintf("r = %.2f", r))
print(subject_cors, n = Inf)
# Read by the manuscript's setup chunk
write_csv(subject_cors |> select(model, parameter, r),
          here("data", "bias_recovery_subject_correlations.csv"))

model_colors <- setNames(eam_palette[c(1, 2, 4)], model_levels)

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
    facet_wrap(~ parameter, nrow = 1, scales = "free", drop = TRUE) +
    # six narrow panels: fewer, smaller tick labels so they do not overlap
    # (the theme() must follow clean_plot(), which would otherwise reset it)
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3) +
    labs(title = model_name,
         x = if (x_title) "True value" else NULL,
         y = if (y_title) "Recovered estimate" else NULL) +
    clean_plot() +
    theme(axis.text = element_text(size = 6))
}

p_subject <- recovery_row("EZ-Diffusion (4par)") /
  recovery_row("EZ-Diffusion (3par)", y_title = TRUE) /
  recovery_row("Full DDM", x_title = TRUE)

ggsave(here("figures", "22_bias_recovery_subject.pdf"),
       p_subject, width = 6.5, height = 7.5)
