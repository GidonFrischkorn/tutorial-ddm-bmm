#' Simulation Study: Evaluating Contaminant Preprocessing Methods
#'
#' Monte Carlo comparison of the three `ezdm_summary_stats()` methods (simple
#' moments, robust median/IQR statistics, and the mixture method — which is the
#' robust-EZ estimator of Wagenmakers, van der Maas, Dolan & Grasman, 2008)
#' under contaminant processes that differ in whether the contaminants overlap
#' the core RT distribution. Structured after ADEMP (Morris et al., 2019) and
#' run with SimDesign (Chalmers & Adkins, 2020).
#'
#' Estimand: the generating DDM parameters (drift, bound, ndt) of every subject
#' x condition cell, recovered by the closed-form EZ equations (Wagenmakers et
#' al., 2007; the same equations bmm's `ezdm()` uses forward), plus two
#' practical estimands (the easy - hard drift effect; subject-level correlations
#' between estimates and truth). The population moments of the core process and
#' an oracle (the plain moments after perfect contaminant exclusion) are carried
#' as diagnostics.
#'
#' Design: 3 contaminant processes x 3 rates (+ shared clean baseline)
#'         x 2 trial counts = 20 rows x 200 replications, 50 subjects each.
#' Contaminant processes (each replaces round(n * p) trials per cell):
#'   fast_guess  anticipations: RT ~ U(0.1 s, ndt_i), response at chance
#'               (Swensson, 1972; Ratcliff, 2008) — outside the core by
#'               construction, i.e. the case the mixture model assumes
#'   delay       delayed start-up: U(0, 2 s) added to a genuine RT, response
#'               kept (Ratcliff & Tuerlinckx, 2002; Ratcliff, 2008;
#'               Wagenmakers et al., 2008) — partly overlaps the core
#'   guess       informationless response: the trial is redrawn from the DDM
#'               with drift = 0 (Vandekerckhove, Tuerlinckx & Lee, 2008;
#'               Chavez De La Pena, Shin & Vandekerckhove, 2026) — fully
#'               overlaps the core, accuracy at chance
#' Methods: simple, robust, mixture (ex-Gaussian / lognormal / inverse-Gaussian
#'          core), oracle; simple, robust and mixture (ex-Gaussian) also after
#'          hard pre-filtering (RT < 150 ms or RT > 5 s). The mixture EM is run
#'          with maxit = 500, under which it converges in > 99.8% of cells; with
#'          bmm's default (maxit = 100) the EM hits the iteration cap in roughly
#'          4-7% of 200-trial cells (18% for the lognormal core) and bmm silently
#'          substitutes robust moments (contaminant_prop = NA).
#' Runtime: ~3.5 h on 15 cores (the ex-Gaussian EM dominates; its cost is bound
#'          by the number of EM iterations, not by the trial count, and per-process
#'          speed roughly halves under full parallel load).

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, tidyverse, SimDesign)
source(here("functions", "clean_plot.R"))

###############################################################################!
# 1) Simulation Design ---------------------------------------------------------
###############################################################################!

# The clean baseline (p_contam = 0) does not depend on the process, so it is
# run once per trial count under process = "none".
Design <- createDesign(
  process  = c("none", "fast_guess", "delay", "guess"),
  p_contam = c(0, 0.02, 0.05, 0.10),
  n_trials = c(50, 200),
  subset   = (process == "none") == (p_contam == 0)
)

# Fixed objects shared across all conditions (not varied)
fo <- list(
  n_subjects = 50,

  # Hyperparameters: each replication draws population means from these
  hyper_drift_easy = c(mean = 2.5, sd = 0.5),
  hyper_drift_hard = c(mean = 1.0, sd = 0.3),
  hyper_bound      = c(mean = 1.5, sd = 0.2),
  hyper_ndt        = c(mean = 0.3, sd = 0.05),

  # Between-subject variability (random effect SDs)
  tau_drift = 0.5,
  tau_bound = 0.15,
  tau_ndt   = 0.03,

  # Contaminant processes
  fast_guess_min = 0.1,     # fast guesses ~ U(0.1 s, ndt_i)
  delay_range    = c(0, 2), # delayed start-up: RT + U(0, 2 s)

  # Hard pre-filtering cutoffs
  hard_lower = 0.15,
  hard_upper = 5.0,

  # Core distributions of the mixture method, and the EM iteration cap. bmm's
  # default (maxit = 100) is hit in a few percent of cells at tol = 1e-6, in
  # which case bmm falls back to robust moments with a warning; with 500 the
  # EM converges in > 99.8% of cells (the remaining fallbacks are counted).
  distributions = c("exgaussian", "lognormal", "invgaussian"),
  em_maxit      = 500
)

###############################################################################!
# 2) EZ Helper Functions -------------------------------------------------------
###############################################################################!

# Forward EZ equations (Wagenmakers et al., 2007, Eqs. 1-4; identical to
# bmm's ezdm_3par_functions.stan): the population proportion of upper
# responses, mean RT and RT variance of a pure Wiener process with s = 1.
# Used for the moment-space diagnostics.
ez_forward <- function(drift, bound, ndt, s = 1) {
  y   <- -drift * bound / s^2
  pc  <- 1 / (1 + exp(y))
  mdt <- (bound / (2 * drift)) * (1 - exp(y)) / (1 + exp(y))
  vrt <- (bound * s^2 / (2 * drift^3)) *
    (2 * y * exp(y) - exp(2 * y) + 1) / (exp(y) + 1)^2
  tibble(pc_pop = pc, mean_rt_pop = ndt + mdt, var_rt_pop = vrt)
}

# Closed-form EZ inversion (Wagenmakers et al., 2007, Eqs. 5-9), s = 1, with
# the standard edge correction for Pc in {0, .5, 1}. Vectorised over cells.
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
# 3) Generate Function ---------------------------------------------------------
###############################################################################!

Generate <- function(condition, fixed_objects) {
  fo   <- fixed_objects
  n_tr <- condition$n_trials

  # Draw population-level parameters for this replication
  pop_drift_easy <- rnorm(1, fo$hyper_drift_easy["mean"],
                          fo$hyper_drift_easy["sd"])
  pop_drift_hard <- max(0.2, rnorm(1, fo$hyper_drift_hard["mean"],
                                   fo$hyper_drift_hard["sd"]))
  pop_bound <- max(0.5, rnorm(1, fo$hyper_bound["mean"],
                              fo$hyper_bound["sd"]))
  pop_ndt   <- max(0.1, rnorm(1, fo$hyper_ndt["mean"],
                              fo$hyper_ndt["sd"]))

  # Subject-level parameters: one random effect per subject, shared across
  # conditions, so the true easy - hard drift effect equals the population
  # difference for every subject. Drift is floored at 0.1 so that every cell
  # is above chance in accuracy coding.
  subj_params <- expand_grid(
    id = 1:fo$n_subjects,
    condition = c("easy", "hard")
  ) |>
    group_by(id) |>
    mutate(
      drift = pmax(0.1, if_else(condition == "easy", pop_drift_easy,
                                pop_drift_hard) + rnorm(1, 0, fo$tau_drift)),
      bound = pop_bound * exp(rnorm(1, 0, fo$tau_bound)),
      ndt   = max(0.1, pop_ndt + rnorm(1, 0, fo$tau_ndt))
    ) |>
    ungroup()

  # Generate clean trial-level data using rddm()
  data_clean <- subj_params |>
    rowwise() |>
    reframe(
      id = id, condition = condition,
      trial = 1:n_tr,
      rddm(n_tr, drift = drift, bound = bound, ndt = ndt, zr = 0.5)
    )

  # Inject contaminants: round(n * p) trials per subject x condition cell
  k <- round(n_tr * condition$p_contam)
  data <- data_clean |>
    left_join(subj_params |> select(id, condition, bound, ndt),
              by = c("id", "condition")) |>
    group_by(id, condition) |>
    mutate(is_contaminant = sample(n()) <= k) |>
    ungroup()

  idx <- which(data$is_contaminant)
  if (length(idx) > 0) {
    if (condition$process == "fast_guess") {
      data$rt[idx]       <- runif(length(idx), fo$fast_guess_min, data$ndt[idx])
      data$response[idx] <- rbinom(length(idx), 1, 0.5)
    } else if (condition$process == "delay") {
      data$rt[idx] <- data$rt[idx] +
        runif(length(idx), fo$delay_range[1], fo$delay_range[2])
    } else if (condition$process == "guess") {
      redraw <- rddm(length(idx), drift = 0, bound = data$bound[idx],
                     ndt = data$ndt[idx], zr = 0.5)
      data$rt[idx]       <- redraw$rt
      data$response[idx] <- redraw$response
    }
  }

  # Truth carried per subject x condition: generating parameters and the
  # population moments of the core process
  truth <- subj_params |>
    rename(drift_true = drift, bound_true = bound, ndt_true = ndt) |>
    mutate(ez_forward(drift_true, bound_true, ndt_true))

  list(
    data  = data |> select(id, condition, trial, rt, response, is_contaminant),
    truth = truth,
    pop   = c(drift_easy = pop_drift_easy, drift_hard = pop_drift_hard,
              bound = pop_bound, ndt = pop_ndt)
  )
}

###############################################################################!
# 4) Analyse Function ----------------------------------------------------------
###############################################################################!

Analyse <- function(condition, dat, fixed_objects) {
  fo    <- fixed_objects
  d     <- dat$data
  truth <- dat$truth
  pop   <- dat$pop

  # Hard pre-filtered version (fixed cutoffs applied before any method)
  d_pf <- d |> filter(rt >= fo$hard_lower, rt <= fo$hard_upper)

  # Summary statistics per method. EM warnings from the mixture method are
  # suppressed here; non-convergence (bmm falls back to robust moments) is
  # tracked below through a missing contaminant proportion.
  sumstats <- function(data, method, label, dist = "exgaussian",
                       maxit = fo$em_maxit) {
    data |>
      summarize(
        suppressWarnings(
          ezdm_summary_stats(rt, response, version = "3par",
                             method = method, distribution = dist,
                             maxit = maxit)
        ),
        .by = c(id, condition)
      ) |>
      mutate(method = label)
  }

  ss <- bind_rows(
    sumstats(d, "simple", "simple"),
    sumstats(d, "robust", "robust"),
    map_dfr(fo$distributions,
            \(dist) sumstats(d, "mixture", paste0("mixture_", dist), dist)),
    sumstats(d |> filter(!is_contaminant), "simple", "oracle"),
    sumstats(d_pf, "simple", "simple_pf"),
    sumstats(d_pf, "robust", "robust_pf"),
    sumstats(d_pf, "mixture", "mixture_exgaussian_pf", "exgaussian")
  )

  # Sanity screen, applied to every method alike: a fitted mean outside the
  # cell's observed RT range or a fitted variance larger than the squared range
  # is a degenerate solution (in pilots the ex-Gaussian EM ran away in about one
  # cell in 10^5, producing e.g. a mean RT of -780 s). Such cells are set to NA
  # and counted, so that a single run-away does not dominate a mean bias.
  rt_range <- d |>
    summarize(rt_min = min(rt), rt_max = max(rt), .by = c(id, condition))

  ss <- ss |>
    left_join(rt_range, by = c("id", "condition")) |>
    mutate(
      degenerate = !is.na(mean_rt) &
        (mean_rt < rt_min | mean_rt > rt_max | var_rt > (rt_max - rt_min)^2),
      mean_rt = if_else(degenerate, NA_real_, mean_rt),
      var_rt  = if_else(degenerate, NA_real_, var_rt)
    ) |>
    select(-rt_min, -rt_max)

  # Closed-form EZ inversion of every cell, joined with the truth and with the
  # oracle's estimates for the same cell
  est <- ss |>
    mutate(ez_invert(mean_rt, var_rt, n_upper, n_trials),
           pc = n_upper / n_trials) |>
    left_join(truth, by = c("id", "condition"))

  oracle_ref <- est |>
    filter(method == "oracle") |>
    select(id, condition,
           drift_orc = drift, bound_orc = bound, ndt_orc = ndt,
           mean_rt_orc = mean_rt, var_rt_orc = var_rt)

  per_cell <- est |>
    left_join(oracle_ref, by = c("id", "condition")) |>
    mutate(
      # errors against the generating parameters (primary estimand)
      e_drift = drift - drift_true,
      e_bound = bound - bound_true,
      e_ndt   = ndt - ndt_true,
      # errors against the population moments of the core (diagnostic)
      e_mean_rt = mean_rt - mean_rt_pop,
      e_var_rt  = var_rt - var_rt_pop,
      e_pc      = pc - pc_pop,
      # errors against the oracle (perfect exclusion) for the same cell
      o_drift   = drift - drift_orc,
      o_bound   = bound - bound_orc,
      o_ndt     = ndt - ndt_orc,
      o_mean_rt = mean_rt - mean_rt_orc,
      o_var_rt  = var_rt - var_rt_orc
    )

  # Per method x condition: mean bias and mean squared error across subjects,
  # subject-level correlations with the truth, mixture diagnostics
  by_cond <- per_cell |>
    summarize(
      bias_drift = mean(e_drift, na.rm = TRUE),
      mse_drift  = mean(e_drift^2, na.rm = TRUE),
      bias_bound = mean(e_bound, na.rm = TRUE),
      mse_bound  = mean(e_bound^2, na.rm = TRUE),
      bias_ndt   = mean(e_ndt, na.rm = TRUE),
      mse_ndt    = mean(e_ndt^2, na.rm = TRUE),
      bias_mean_rt = mean(e_mean_rt, na.rm = TRUE),
      bias_var_rt  = mean(e_var_rt, na.rm = TRUE),
      bias_pc      = mean(e_pc, na.rm = TRUE),
      orc_drift   = mean(o_drift, na.rm = TRUE),
      orc_bound   = mean(o_bound, na.rm = TRUE),
      orc_ndt     = mean(o_ndt, na.rm = TRUE),
      orc_mean_rt = mean(o_mean_rt, na.rm = TRUE),
      orc_var_rt  = mean(o_var_rt, na.rm = TRUE),
      cor_drift = suppressWarnings(cor(drift, drift_true, use = "complete.obs")),
      cor_bound = suppressWarnings(cor(bound, bound_true, use = "complete.obs")),
      cor_ndt   = suppressWarnings(cor(ndt, ndt_true, use = "complete.obs")),
      pi_c       = mean(contaminant_prop, na.rm = TRUE),
      fallback   = mean(is.na(contaminant_prop)),
      degenerate = mean(degenerate),
      na_rate    = mean(is.na(drift)),
      .by = c(method, condition)
    )

  # Easy - hard drift effect: the mean over subjects of the estimated difference
  # against the mean of the subjects' generating differences (these differ from
  # the population difference only where a subject's hard-condition drift was
  # floored at 0.1)
  effect <- est |>
    select(method, id, condition, drift, drift_true) |>
    pivot_wider(names_from = condition, values_from = c(drift, drift_true)) |>
    summarize(
      bias_effect = mean(drift_easy - drift_hard, na.rm = TRUE) -
        mean(drift_true_easy - drift_true_hard, na.rm = TRUE),
      .by = method
    ) |>
    mutate(mse_effect = bias_effect^2, condition = "all")

  # Return as a named vector: <measure>__<method>__<condition>. The set of
  # names must be identical in every replication, so the mixture diagnostics
  # are only returned for the mixture methods and the effect measures only
  # once per method.
  out <- bind_rows(
    by_cond |>
      pivot_longer(-c(method, condition), names_to = "measure") |>
      filter(!(measure %in% c("pi_c", "fallback") &
                 !str_starts(method, "mixture"))),
    effect |>
      pivot_longer(c(bias_effect, mse_effect), names_to = "measure")
  ) |>
    mutate(key = paste(measure, method, condition, sep = "__"))

  setNames(out$value, out$key)
}

###############################################################################!
# 5) Summarise Function --------------------------------------------------------
###############################################################################!

Summarise <- function(condition, results, fixed_objects) {
  # results: matrix with one row per replication and one named column per
  # statistic. Bias-type measures are averaged; MSE columns become RMSE; the
  # Monte Carlo SE is reported for every bias measure.
  m <- colMeans(results, na.rm = TRUE)
  is_mse <- str_starts(names(m), "mse_")
  m[is_mse] <- sqrt(m[is_mse])
  names(m)[is_mse] <- str_replace(names(m)[is_mse], "^mse_", "rmse_")

  is_bias <- str_starts(names(results), "bias_")
  mcse <- apply(results[, is_bias, drop = FALSE], 2, sd, na.rm = TRUE) /
    sqrt(colSums(!is.na(results[, is_bias, drop = FALSE])))
  names(mcse) <- str_replace(names(mcse), "^bias_", "mcse_")

  c(m, mcse)
}

###############################################################################!
# 6) Run Simulation ------------------------------------------------------------
###############################################################################!

n_reps  <- 200
n_cores <- min(parallel::detectCores() - 1, 15)

cat("=== Simulation Design ===\n")
cat("Conditions:", nrow(Design), "\n")
cat("Replications per condition:", n_reps, "\n")
cat("Total datasets:", nrow(Design) * n_reps, "\n")
cat("Parallel cores:", n_cores, "\n\n")

# The same seed for every design row gives common random numbers: within a
# trial count, the clean data are identical across processes and rates, so
# differences between rows are due to the contaminants, not to sampling.
sim_results <- runSimulation(
  design        = Design,
  replications  = n_reps,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = fo,
  packages      = c("bmm", "tidyverse"),
  parallel      = TRUE,
  ncores        = n_cores,
  save          = TRUE,
  save_results  = TRUE,
  filename      = "output/simulation_preprocessing_methods",
  seed          = rep(2026, nrow(Design)),
  progress      = TRUE,
  # a missing correlation or contaminant proportion in one replication is a
  # legitimate result, not a reason to redraw the data
  control       = list(allow_na = TRUE, allow_nan = TRUE)
)

###############################################################################!
# 7) Tidy Results --------------------------------------------------------------
###############################################################################!

method_levels <- c("simple", "robust", "mixture_exgaussian",
                   "mixture_lognormal", "mixture_invgaussian", "oracle",
                   "simple_pf", "robust_pf", "mixture_exgaussian_pf")
method_labels <- c("Simple", "Robust", "Mixture (ex-Gaussian)",
                   "Mixture (lognormal)", "Mixture (inv-Gaussian)",
                   "Oracle (perfect exclusion)",
                   "Simple, pre-filtered", "Robust, pre-filtered",
                   "Mixture (ex-Gaussian), pre-filtered")
process_labels <- c(fast_guess = "Fast guesses (below ndt)",
                    delay = "Delayed start-ups (+ U(0, 2 s))",
                    guess = "Informationless responses (drift 0)")

results_long <- sim_results |>
  as_tibble() |>
  select(process, p_contam, n_trials, matches("__")) |>
  pivot_longer(matches("__"), names_to = "key", values_to = "value") |>
  separate_wider_delim(key, "__", names = c("measure", "method", "condition"))

# The clean baseline row is shared by all processes; replicate it per process
# so that every process has a p = 0 point.
baseline <- results_long |>
  filter(process == "none") |>
  select(-process) |>
  expand_grid(process = names(process_labels))

results_long <- bind_rows(results_long |> filter(process != "none"), baseline) |>
  mutate(
    method = factor(method, levels = method_levels, labels = method_labels),
    process = factor(process, levels = names(process_labels),
                     labels = process_labels),
    condition = factor(condition, levels = c("easy", "hard", "all"),
                       labels = c("Easy", "Hard", "All"))
  ) |>
  arrange(n_trials, process, p_contam, method, condition, measure)

write_csv(results_long, here("data", "simulation_preprocessing_summary.csv"))

###############################################################################!
# 8) Figures -------------------------------------------------------------------
###############################################################################!

param_labels <- c(drift = "Drift", bound = "Bound", ndt = "ndt")

# Facet strips are narrow: break each process label before its parenthesis
wrap_process <- \(x) str_replace(x, " \\(", "\n(")

## 8.1) Figure 1: signed bias in drift, bound, ndt by process (n = 200) --------

main_methods <- c("Simple", "Robust", "Mixture (ex-Gaussian)",
                  "Oracle (perfect exclusion)")
main_colors <- c(
  "Simple" = "grey50",
  "Robust" = eam_palette[2],
  "Mixture (ex-Gaussian)" = eam_palette[1],
  "Oracle (perfect exclusion)" = "black"
)

fig1_data <- results_long |>
  filter(n_trials == 200, method %in% main_methods,
         measure %in% c("bias_drift", "bias_bound", "bias_ndt")) |>
  mutate(parameter = factor(str_remove(measure, "bias_"),
                            levels = names(param_labels),
                            labels = param_labels))

p_fig1 <- fig1_data |>
  ggplot(aes(x = p_contam, y = value, color = method, linetype = condition,
             group = interaction(method, condition))) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_grid(parameter ~ process, scales = "free_y",
             labeller = labeller(process = wrap_process)) +
  scale_color_manual(values = main_colors) +
  scale_linetype_manual(values = c("Easy" = "solid", "Hard" = "dashed")) +
  scale_x_continuous(labels = scales::percent, breaks = c(0, .02, .05, .10)) +
  labs(x = "Contamination rate", y = "Bias (estimate − generating value)",
       color = "Method", linetype = "Condition") +
  clean_plot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 8), legend.spacing.y = unit(0, "pt"))

ggsave(here("figures", "02_method_comparison.pdf"),
       p_fig1, width = 6.5, height = 6, device = "pdf")

## 8.2) Figure 2 (supplementary): mixture core distributions (n = 200) --------

core_methods <- c("Mixture (ex-Gaussian)", "Mixture (lognormal)",
                  "Mixture (inv-Gaussian)")
core_colors <- c(
  "Mixture (ex-Gaussian)" = eam_palette[1],
  "Mixture (lognormal)" = eam_palette[3],
  "Mixture (inv-Gaussian)" = eam_palette[4]
)

fig2_data <- results_long |>
  filter(n_trials == 200, method %in% core_methods,
         measure %in% c("bias_drift", "bias_bound", "bias_ndt")) |>
  mutate(parameter = factor(str_remove(measure, "bias_"),
                            levels = names(param_labels),
                            labels = param_labels))

p_fig2 <- fig2_data |>
  ggplot(aes(x = p_contam, y = value, color = method, linetype = condition,
             group = interaction(method, condition))) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_grid(parameter ~ process, scales = "free_y",
             labeller = labeller(process = wrap_process)) +
  scale_color_manual(values = core_colors) +
  scale_linetype_manual(values = c("Easy" = "solid", "Hard" = "dashed")) +
  scale_x_continuous(labels = scales::percent, breaks = c(0, .02, .05, .10)) +
  labs(x = "Contamination rate", y = "Bias (estimate − generating value)",
       color = "Core distribution", linetype = "Condition") +
  clean_plot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 8), legend.spacing.y = unit(0, "pt"))

ggsave(here("figures", "02_mixture_core_distributions.pdf"),
       p_fig2, width = 6.5, height = 6, device = "pdf")

###############################################################################!
# 9) Console Summary -----------------------------------------------------------
###############################################################################!

options(width = 200)

cat("\n=== Clean data (p = 0): bias in EZ parameters and correlations ===\n\n")
results_long |>
  filter(p_contam == 0, process == process_labels[1], condition != "All",
         measure %in% c("bias_drift", "bias_bound", "bias_ndt",
                        "cor_drift", "cor_bound", "cor_ndt")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(n_trials, method, condition, starts_with("bias_"), starts_with("cor_")) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== 5% contamination, n = 200: bias in EZ parameters by process ===\n\n")
results_long |>
  filter(p_contam == 0.05, n_trials == 200, condition != "All",
         measure %in% c("bias_drift", "bias_bound", "bias_ndt", "pi_c")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(process, method, condition, bias_drift, bias_bound, bias_ndt, pi_c) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== 5% contamination, n = 200: practical estimands ===\n\n")
results_long |>
  filter(p_contam == 0.05, n_trials == 200,
         measure %in% c("bias_effect", "rmse_effect", "cor_drift",
                        "cor_bound", "cor_ndt")) |>
  pivot_wider(names_from = c(measure, condition), values_from = value) |>
  select(process, method, bias_effect_All, rmse_effect_All,
         cor_drift_Easy, cor_drift_Hard, cor_bound_Easy, cor_ndt_Easy) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== Pre-filtering: bias vs the unfiltered method (5%, n = 200) ===\n\n")
results_long |>
  filter(p_contam == 0.05, n_trials == 200, condition != "All",
         measure %in% c("bias_drift", "bias_bound", "bias_ndt"),
         str_detect(method, "Simple|Robust|Mixture \\(ex-Gaussian\\)")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(process, method, condition, bias_drift, bias_bound, bias_ndt) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== EM non-convergence at maxit = 500 (fallback) and degenerate fits ===\n\n")
results_long |>
  filter(measure %in% c("fallback", "degenerate"), condition != "All") |>
  summarize(rate = signif(mean(value), 3), .by = c(measure, n_trials, method)) |>
  pivot_wider(names_from = measure, values_from = rate) |>
  print(n = Inf)
