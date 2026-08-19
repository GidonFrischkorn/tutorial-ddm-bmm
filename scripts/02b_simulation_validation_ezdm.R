#' Validation: Hierarchical ezdm() Fits of Simulated Contaminated Data
#'
#' The Monte Carlo study in scripts/02_simulation_contaminant_methods.R
#' evaluates the preprocessing methods with the closed-form EZ equations. This
#' script checks that its conclusions carry into the hierarchical Bayesian
#' `ezdm()` fits the tutorial actually uses. It regenerates one dataset per
#' contaminant process (5% contamination, 200 trials, 50 subjects; the same
#' generating mechanism as script 02, with population values fixed at the
#' hyperparameter means) plus one clean dataset, computes the summary statistics
#' with the simple, robust and mixture (ex-Gaussian) methods, fits `ezdm()` to
#' each, and compares posterior means with the generating values and with the
#' closed-form estimates.
#'
#' 4 datasets x 3 methods = 12 fits (~20 min with 5 chains on 5 cores).

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, brms, tidyverse)
source(here("functions", "extract_fixef.R"))

# -- Fitting settings (as in the tutorials) -----------------------------------
chains <- 5
cores  <- 5
warmup <- 1000
iter   <- 2000
adapt_delta   <- 0.90
max_treedepth <- 12

# -- Closed-form EZ inversion (Wagenmakers et al., 2007), as in script 02 -----
ez_invert <- function(mean_rt, var_rt, n_upper, n_trials, s = 1) {
  pc <- n_upper / n_trials
  pc <- ifelse(pc == 1,   1 - 1 / (2 * n_trials), pc)
  pc <- ifelse(pc == 0,   1 / (2 * n_trials), pc)
  pc <- ifelse(pc == 0.5, 0.5 + 1 / (2 * n_trials), pc)
  L  <- qlogis(pc)
  x  <- L * (L * pc^2 - L * pc + pc - 0.5) / var_rt
  drift <- sign(pc - 0.5) * s * x^(1 / 4)
  bound <- s^2 * L / drift
  y     <- -drift * bound / s^2
  mdt   <- (bound / (2 * drift)) * (1 - exp(y)) / (1 + exp(y))
  tibble(drift = drift, bound = bound, ndt = mean_rt - mdt)
}

###############################################################################!
# 1) Simulate One Dataset per Contaminant Process ------------------------------
###############################################################################!

set.seed(2026)
n_subjects <- 50
n_trials   <- 200
p_contam   <- 0.05

# Subject-level parameters (mirrors script 02: one random effect per subject
# shared across conditions; population values at the hyperparameter means)
subj_params <- expand_grid(id = 1:n_subjects, condition = c("easy", "hard")) |>
  group_by(id) |>
  mutate(
    drift = pmax(0.1, if_else(condition == "easy", 2.5, 1.0) + rnorm(1, 0, 0.5)),
    bound = 1.5 * exp(rnorm(1, 0, 0.15)),
    ndt   = max(0.1, 0.3 + rnorm(1, 0, 0.03))
  ) |>
  ungroup()

data_clean <- subj_params |>
  rowwise() |>
  reframe(
    id = id, condition = condition, trial = 1:n_trials,
    rddm(n_trials, drift = drift, bound = bound, ndt = ndt, zr = 0.5)
  ) |>
  left_join(subj_params |> select(id, condition, bound, ndt),
            by = c("id", "condition")) |>
  group_by(id, condition) |>
  mutate(is_contaminant = sample(n()) <= round(n() * p_contam)) |>
  ungroup()

# The same trial positions are contaminated under each process
idx <- which(data_clean$is_contaminant)
data_fast <- data_clean
data_fast$rt[idx]       <- runif(length(idx), 0.1, data_clean$ndt[idx])
data_fast$response[idx] <- rbinom(length(idx), 1, 0.5)

data_delay <- data_clean
data_delay$rt[idx] <- data_clean$rt[idx] + runif(length(idx), 0, 2)

data_guess <- data_clean
redraw <- rddm(length(idx), drift = 0, bound = data_clean$bound[idx],
               ndt = data_clean$ndt[idx], zr = 0.5)
data_guess$rt[idx]       <- redraw$rt
data_guess$response[idx] <- redraw$response

datasets <- list(none = data_clean, fast_guess = data_fast,
                 delay = data_delay, guess = data_guess)

###############################################################################!
# 2) Summary Statistics and ezdm() Fits ----------------------------------------
###############################################################################!

model_ezdm <- ezdm(mean_rt = "mean_rt", var_rt = "var_rt",
                   n_upper = "n_upper", n_trials = "n_trials",
                   version = "3par")

# Condition-specific drift and bound (as in Tutorial 1), so that a spurious
# bound difference induced by contamination can show up; ndt is shared
formula_ezdm <- bmf(
  drift ~ 0 + condition + (0 + condition | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

methods <- c("simple", "robust", "mixture")

sumstats_all <- expand_grid(process = names(datasets), method = methods) |>
  mutate(data = map2(process, method, \(pr, me) {
    datasets[[pr]] |>
      summarize(
        suppressWarnings(
          ezdm_summary_stats(rt, response, version = "3par", method = me,
                             distribution = "exgaussian", maxit = 500)
        ),
        .by = c(id, condition)
      )
  }))

fits <- sumstats_all |>
  mutate(fit = pmap(list(process, method, data), \(pr, me, d) {
    bmm(
      formula = formula_ezdm,
      data    = d,
      model   = model_ezdm,
      chains  = chains, cores = cores, warmup = warmup, iter = iter,
      backend = "cmdstanr", sample_prior = "yes",
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      file    = here("output", paste0("02b_ezdm_", pr, "_", me))
    )
  }))

###############################################################################!
# 3) Compare Posterior Means with Truth and Closed-Form Estimates --------------
###############################################################################!

# Truth: the mean of the subject-level generating values in this dataset
truth <- bind_rows(
  subj_params |>
    summarize(truth = mean(drift), .by = condition) |>
    mutate(parameter = "drift", level = condition),
  subj_params |>
    summarize(truth = mean(bound), .by = condition) |>
    mutate(parameter = "bound", level = condition),
  subj_params |>
    summarize(truth = mean(ndt)) |>
    mutate(parameter = "ndt", level = "overall")
) |>
  select(parameter, level, truth)

# Closed-form estimates: EZ inversion of each cell, averaged over subjects
closed_form <- sumstats_all |>
  mutate(cf = map(data, \(d) {
    est <- d |> mutate(ez_invert(mean_rt, var_rt, n_upper, n_trials))
    bind_rows(
      est |> summarize(closed_form = mean(drift, na.rm = TRUE), .by = condition) |>
        mutate(parameter = "drift", level = condition),
      est |> summarize(closed_form = mean(bound, na.rm = TRUE), .by = condition) |>
        mutate(parameter = "bound", level = condition),
      est |> summarize(closed_form = mean(ndt, na.rm = TRUE)) |>
        mutate(parameter = "ndt", level = "overall")
    ) |>
      select(parameter, level, closed_form)
  })) |>
  select(process, method, cf) |>
  unnest(cf)

# Hierarchical fits: fixed effects on the parameter scale
posterior <- fits |>
  mutate(fe = map(fit, \(f) {
    extract_fixef(f, "ezdm", c("drift", "bound", "ndt"),
                  list(bound = exp, ndt = exp)) |>
      clean_levels() |>
      mutate(level = if_else(level == "", "overall", level)) |>
      select(parameter, level, posterior_mean = estimate, lower, upper)
  })) |>
  select(process, method, fe) |>
  unnest(fe)

validation <- posterior |>
  left_join(closed_form, by = c("process", "method", "parameter", "level")) |>
  left_join(truth, by = c("parameter", "level")) |>
  mutate(
    bias_posterior   = posterior_mean - truth,
    bias_closed_form = closed_form - truth,
    process = factor(process, levels = names(datasets)),
    method  = factor(method, levels = methods),
    parameter = factor(parameter, levels = c("drift", "bound", "ndt"))
  ) |>
  arrange(process, method, parameter, level) |>
  select(process, method, parameter, level, truth, closed_form,
         posterior_mean, lower, upper, bias_closed_form, bias_posterior)

write_csv(validation, here("data", "simulation_validation_ezdm.csv"))

# Drift effect (easy - hard) from the posterior draws vs truth
effect <- fits |>
  mutate(eff = map(fit, \(f) {
    draws <- as_draws_df(f)
    e <- draws$b_drift_conditioneasy - draws$b_drift_conditionhard
    tibble(effect_mean = mean(e), lower = quantile(e, .025), upper = quantile(e, .975))
  })) |>
  select(process, method, eff) |>
  unnest(eff) |>
  mutate(effect_true = with(subj_params, mean(drift[condition == "easy"]) -
                              mean(drift[condition == "hard"])),
         bias = effect_mean - effect_true)

options(width = 200)
cat("\n=== Hierarchical ezdm() vs closed-form: bias against generating values ===\n\n")
validation |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== Easy - hard drift effect (posterior mean and 95% CrI) ===\n\n")
effect |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

# Convergence check
cat("\n=== Max Rhat per fit ===\n\n")
fits |>
  mutate(max_rhat = map_dbl(fit, \(f) max(brms::rhat(f), na.rm = TRUE))) |>
  select(process, method, max_rhat) |>
  print(n = Inf)
