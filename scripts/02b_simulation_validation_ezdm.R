#' Validation: Hierarchical ezdm() Fits of Simulated Contaminated Data
#'
#' The Monte Carlo study in scripts/02_simulation_contaminant_methods.R
#' evaluates the preprocessing methods with the closed-form EZ equations. This
#' script checks that its conclusions carry into the hierarchical Bayesian
#' `ezdm()` fits the tutorial actually uses. It regenerates one dataset per
#' contaminant process of the v7 design — leading-edge anticipations (depth
#' x_i ~ U(0.1, 0.5), RT ~ U((1 - x_i) ndt_i, (1 + x_i) ndt_i), response at
#' chance), delayed start-ups, and zero-drift guesses — at 5% contamination,
#' 200 trials, 50 subjects (the same generating mechanism as script 02, with
#' population values fixed at the hyperparameter means) plus one clean
#' dataset, computes the summary statistics with the simple, robust and
#' mixture (ex-Gaussian) methods plus the accuracy-adjusted mixture
#' (mixture statistics with adjust_ezdm_accuracy() applied to the response
#' counts), fits `ezdm()` to each, and compares posterior means with the
#' generating values and with the closed-form estimates.
#'
#' In addition to the pure processes, one dataset of the mixed process — the
#' primary row of script 02: each subject draws a total rate
#' p_i ~ Beta(mean .05, concentration 10) and a process mix
#' w_i ~ Dirichlet(1, 1, 1) — validates the bridge for the row that carries
#' the recommendation. The mixed row is the only one in which hierarchical
#' partial pooling can interact with contamination, because only there do
#' subjects differ in how much and how they are contaminated: shrinkage pulls
#' the heavily contaminated minority toward a group mean the near-clean
#' majority dominates. Section 4 therefore compares the closed-form and the
#' hierarchical subject-level drift errors against the subjects' true rates.
#' Ex-ante prediction (recorded 2026-08-24, before the mixed fits were run):
#' shrinkage compresses subject-level errors but cannot recreate the
#' information contamination destroyed, so the error-rate correlation should
#' persist at similar or modestly reduced magnitude in the hierarchical fits.
#' If it were substantially repaired instead, the manuscript's design-side
#' claim ("no preprocessing choice repaired the differential bias") would be
#' overclaimed for `ezdm` users and would need softening.
#'
#' 5 datasets x 4 methods = 20 fits (~30 min with 5 chains on 5 cores).

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

# Leading-edge anticipations: each subject draws a depth x_i ~ U(0.1, 0.5)
# (shared across conditions); RT ~ U((1 - x_i) ndt_i, (1 + x_i) ndt_i),
# response at chance (as in script 02, v7)
subj_x <- tibble(id = 1:n_subjects, x_i = runif(n_subjects, 0.1, 0.5))
data_led <- data_clean |> left_join(subj_x, by = "id")
data_led$rt[idx]       <- runif(length(idx),
                                (1 - data_led$x_i[idx]) * data_led$ndt[idx],
                                (1 + data_led$x_i[idx]) * data_led$ndt[idx])
data_led$response[idx] <- rbinom(length(idx), 1, 0.5)
data_led <- data_led |> select(-x_i)

data_delay <- data_clean
data_delay$rt[idx] <- data_clean$rt[idx] + runif(length(idx), 0, 2)

data_guess <- data_clean
redraw <- rddm(length(idx), drift = 0, bound = data_clean$bound[idx],
               ndt = data_clean$ndt[idx], zr = 0.5)
data_guess$rt[idx]       <- redraw$rt
data_guess$response[idx] <- redraw$response

# Mixed process (as in script 02): subject-specific rates and process mixes,
# so this is the one dataset in which contamination varies between subjects.
# Generated under a local seed with the surrounding RNG state restored
# afterwards, so the stochastic adjust_ezdm_accuracy() draws for the
# pure-process datasets below consume exactly the stream positions they
# consumed before this section existed — otherwise the cached pure-process
# mixture_adjust fits would silently mismatch their regenerated data.
seed_before <- get(".Random.seed", envir = .GlobalEnv)
set.seed(2027)
types <- c("leading_edge", "delay", "guess")
w <- matrix(rexp(n_subjects * 3), ncol = 3, dimnames = list(NULL, types))
profile_mixed <- tibble(
  id  = 1:n_subjects,
  p_i = rbeta(n_subjects, p_contam * 10, (1 - p_contam) * 10),
  x_i = runif(n_subjects, 0.1, 0.5)
) |>
  bind_cols(as_tibble(w / rowSums(w)))

data_mixed <- data_clean |>
  select(-is_contaminant) |>
  left_join(profile_mixed, by = "id") |>
  group_by(id, condition) |>
  mutate(is_contaminant = sample(n()) <= round(n() * first(p_i))) |>
  ungroup()

idx_m <- which(data_mixed$is_contaminant)
contam_type <- vapply(idx_m, function(i) {
  sample(types, 1, prob = c(data_mixed$leading_edge[i], data_mixed$delay[i],
                            data_mixed$guess[i]))
}, character(1))

i_le <- idx_m[contam_type == "leading_edge"]
data_mixed$rt[i_le]       <- runif(length(i_le),
                                   (1 - data_mixed$x_i[i_le]) * data_mixed$ndt[i_le],
                                   (1 + data_mixed$x_i[i_le]) * data_mixed$ndt[i_le])
data_mixed$response[i_le] <- rbinom(length(i_le), 1, 0.5)
i_dl <- idx_m[contam_type == "delay"]
data_mixed$rt[i_dl] <- data_mixed$rt[i_dl] + runif(length(i_dl), 0, 2)
i_gu <- idx_m[contam_type == "guess"]
redraw_m <- rddm(length(i_gu), drift = 0, bound = data_mixed$bound[i_gu],
                 ndt = data_mixed$ndt[i_gu], zr = 0.5)
data_mixed$rt[i_gu]       <- redraw_m$rt
data_mixed$response[i_gu] <- redraw_m$response

data_mixed <- data_mixed |>
  select(-x_i, -leading_edge, -delay, -guess)
assign(".Random.seed", seed_before, envir = .GlobalEnv)

datasets <- list(none = data_clean, leading_edge = data_led,
                 delay = data_delay, guess = data_guess, mixed = data_mixed)

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

methods <- c("simple", "robust", "mixture", "mixture_adjust")

sumstats_base <- expand_grid(process = names(datasets),
                             method = c("simple", "robust", "mixture")) |>
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

# Accuracy-adjusted mixture: the mixture statistics with the response counts
# corrected for the estimated share of chance-accuracy contaminants
# (adjust_ezdm_accuracy() is stochastic and scalar-only, so it is applied per
# cell; a cell with contaminant_prop NA is left unchanged). Only n_upper and
# n_trials change; the RT moments are the mixture's.
sumstats_adj <- sumstats_base |>
  filter(method == "mixture") |>
  mutate(
    method = "mixture_adjust",
    data = map(data, \(d) {
      adj <- pmap(list(d$n_upper, d$n_trials, d$contaminant_prop),
                  \(nu, nt, pc) adjust_ezdm_accuracy(nu, nt, pc)) |>
        bind_rows()
      d |> mutate(n_upper = adj$n_upper_adj, n_trials = adj$n_trials_adj)
    })
  )

sumstats_all <- bind_rows(sumstats_base, sumstats_adj)

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

###############################################################################!
# 4) Differential Bias in the Mixed Row: Hierarchical vs Closed-Form -----------
###############################################################################!

# The closed-form study's headline for the mixed row is differential bias:
# each subject's drift error tracks their own contamination rate under every
# real method. Partial pooling is the one piece of machinery the closed-form
# analysis does not have, so this is the comparison the mixed fits exist for
# (ex-ante prediction in the header). Both sides are computed on the same
# subjects: the closed-form error from the cell's EZ inversion, the
# hierarchical error from the subject-specific posterior means
# (fixed + random effects via coef()).
subj_truth <- subj_params |> select(id, condition, drift_true = drift)

diff_mixed <- fits |>
  filter(process == "mixed") |>
  mutate(dd = pmap(list(method, data, fit), \(me, d, f) {
    cf <- d |>
      mutate(ez_invert(mean_rt, var_rt, n_upper, n_trials)) |>
      select(id, condition, drift_cf = drift)
    co <- coef(f)$id
    keep <- grep("^drift_", dimnames(co)[[3]], value = TRUE)
    hier <- as_tibble(co[, "Estimate", keep], rownames = "id") |>
      mutate(id = as.integer(id)) |>
      pivot_longer(-id, names_to = "condition", values_to = "drift_hier") |>
      mutate(condition = str_remove(condition, "^drift_condition"))
    cf |>
      left_join(hier, by = c("id", "condition")) |>
      left_join(subj_truth, by = c("id", "condition")) |>
      left_join(profile_mixed |> select(id, p_i), by = "id") |>
      summarize(
        cor_pcon_cf   = cor(drift_cf - drift_true, p_i, use = "complete.obs"),
        cor_pcon_hier = cor(drift_hier - drift_true, p_i, use = "complete.obs"),
        sd_error_cf   = sd(drift_cf - drift_true, na.rm = TRUE),
        sd_error_hier = sd(drift_hier - drift_true, na.rm = TRUE),
        .by = condition
      )
  })) |>
  select(method, dd) |>
  unnest(dd)

write_csv(diff_mixed, here("data", "simulation_validation_mixed_differential.csv"))

cat("\n=== Mixed row: subject drift error vs true contamination rate ===\n")
cat("(cor_pcon_*: correlation of the subject's drift error with their true\n")
cat(" rate p_i; sd_error_*: spread of the subject-level errors, where\n")
cat(" shrinkage should show as compression in the hierarchical column)\n\n")
diff_mixed |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)
