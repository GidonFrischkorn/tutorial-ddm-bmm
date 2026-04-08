#' Simulation Study: Evaluating Contaminant Detection Methods
#'
#' Monte Carlo simulation comparing preprocessing methods for EAM data across
#' varying contamination rates. Uses the SimDesign package for structured
#' simulation with Generate/Analyse/Summarise workflow.
#'
#' Design: 3 (p_fastguess) x 3 (p_slowcont) x 500 repetitions = 4,500 datasets
#' Methods: simple, robust, mixture (exgaussian/lognormal/invgaussian),
#'          each with and without hard pre-filtering (RT < 150ms or RT > 5s)

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, tidyverse, SimDesign)
source(here("functions", "clean_plot.R"))

###############################################################################!
# 1) Simulation Design ---------------------------------------------------------
###############################################################################!

Design <- createDesign(
  p_fastguess = c(0, 0.02, 0.05),
  p_slowcont  = c(0, 0.02, 0.05)
)

# Fixed objects shared across all conditions (not varied)
fo <- list(
  n_subjects    = 50,
  n_trials_cond = 200,

  # Hyperparameters: each rep draws population means from these
  hyper_drift_easy = c(mean = 2.5, sd = 0.5),
  hyper_drift_hard = c(mean = 1.0, sd = 0.3),
  hyper_bound      = c(mean = 1.5, sd = 0.2),
  hyper_ndt        = c(mean = 0.3, sd = 0.05),

  # Between-subject variability (random effect SDs)
  tau_drift = 0.5,
  tau_bound = 0.15,
  tau_ndt   = 0.03,

  # Hard pre-filtering cutoffs
  hard_lower = 0.15,
  hard_upper = 5.0,

  # Mixture distributions to evaluate
  distributions = c("exgaussian", "lognormal", "invgaussian")
)

###############################################################################!
# 2) Generate Function ---------------------------------------------------------
###############################################################################!

Generate <- function(condition, fixed_objects) {
  fo <- fixed_objects

  # Draw population-level parameters for this repetition
  pop_drift_easy <- rnorm(1, fo$hyper_drift_easy["mean"],
                          fo$hyper_drift_easy["sd"])
  pop_drift_hard <- max(0.2, rnorm(1, fo$hyper_drift_hard["mean"],
                                   fo$hyper_drift_hard["sd"]))
  pop_bound <- max(0.5, rnorm(1, fo$hyper_bound["mean"],
                              fo$hyper_bound["sd"]))
  pop_ndt   <- max(0.1, rnorm(1, fo$hyper_ndt["mean"],
                               fo$hyper_ndt["sd"]))

  # Generate subject-level parameters
  subj_params <- expand_grid(
    id = 1:fo$n_subjects,
    condition = c("easy", "hard")
  ) |>
    group_by(id) |>
    mutate(
      drift = if_else(condition == "easy", pop_drift_easy, pop_drift_hard) +
        rnorm(1, 0, fo$tau_drift),
      bound = pop_bound * exp(rnorm(1, 0, fo$tau_bound)),
      ndt   = max(0.1, pop_ndt + rnorm(1, 0, fo$tau_ndt))
    ) |>
    ungroup()

  # Generate trial-level data using rddm()
  data_clean <- subj_params |>
    rowwise() |>
    reframe(
      id = id, condition = condition,
      trial = 1:fo$n_trials_cond,
      rddm(fo$n_trials_cond, drift = drift, bound = bound,
            ndt = ndt, zr = 0.5)
    )

  # Compute clean summary stats (ground truth)
  sumstats_clean <- data_clean |>
    summarize(
      ezdm_summary_stats(rt, response, version = "3par", method = "simple"),
      .by = c(id, condition)
    )

  # Inject contaminants
  p_fast <- condition$p_fastguess
  p_slow <- condition$p_slowcont

  if (p_fast == 0 && p_slow == 0) {
    data_contam <- data_clean |> mutate(is_contaminant = FALSE)
  } else {
    rt_stats <- data_clean |>
      group_by(id) |>
      summarize(mean_rt_subj = mean(rt), sd_rt_subj = sd(rt),
                .groups = "drop")

    ndt_by_subj <- subj_params |> distinct(id, ndt)

    data_contam <- data_clean |>
      left_join(rt_stats, by = "id") |>
      left_join(ndt_by_subj, by = "id") |>
      group_by(id, condition) |>
      mutate(
        n = n(),
        contam_draw = sample(n()),
        is_fast = contam_draw <= round(n * p_fast),
        is_slow = !is_fast & contam_draw > (n - round(n * p_slow))
      ) |>
      ungroup() |>
      mutate(
        rt = if_else(is_fast, runif(n(), 0.1, ndt), rt),
        response = if_else(is_fast,
                           as.numeric(rbinom(n(), 1, 0.5)), response),
        rt = if_else(is_slow,
                     runif(n(), mean_rt_subj + 2 * sd_rt_subj,
                           mean_rt_subj + 5 * sd_rt_subj), rt),
        is_contaminant = is_fast | is_slow
      ) |>
      select(id, condition, trial, rt, response, is_contaminant)
  }

  list(
    data_contam = data_contam,
    sumstats_clean = sumstats_clean,
    pop_params = c(drift_easy = pop_drift_easy,
                   drift_hard = pop_drift_hard,
                   bound = pop_bound, ndt = pop_ndt)
  )
}

###############################################################################!
# 3) Analyse Function ----------------------------------------------------------
###############################################################################!

Analyse <- function(condition, dat, fixed_objects) {
  fo <- fixed_objects
  data_contam    <- dat$data_contam
  sumstats_clean <- dat$sumstats_clean

  # Pre-filtered version
  data_hard <- data_contam |>
    filter(rt >= fo$hard_lower & rt <= fo$hard_upper)

  input_data <- list(none = data_contam, hard = data_hard)

  results <- list()

  # --- Simple and robust baselines ---
  for (pf in names(input_data)) {
    d <- input_data[[pf]]

    ss_simple <- d |>
      summarize(
        ezdm_summary_stats(rt, response, version = "3par", method = "simple"),
        .by = c(id, condition)
      )
    results <- c(results, list(tibble(
      approach = "simple", distribution = "none", prefilter = pf,
      id = ss_simple$id, condition = ss_simple$condition,
      mean_rt = ss_simple$mean_rt, var_rt = ss_simple$var_rt
    )))

    ss_robust <- d |>
      summarize(
        ezdm_summary_stats(rt, response, version = "3par", method = "robust"),
        .by = c(id, condition)
      )
    results <- c(results, list(tibble(
      approach = "robust", distribution = "none", prefilter = pf,
      id = ss_robust$id, condition = ss_robust$condition,
      mean_rt = ss_robust$mean_rt, var_rt = ss_robust$var_rt
    )))
  }

  # --- Mixture method across distributions x prefilter ---
  for (dist in fo$distributions) {
    for (pf in names(input_data)) {
      d <- input_data[[pf]]
      ss <- tryCatch(
        d |>
          summarize(
            ezdm_summary_stats(rt, response, version = "3par",
                               method = "mixture", distribution = dist),
            .by = c(id, condition)
          ),
        error = function(e) NULL
      )
      if (!is.null(ss)) {
        results <- c(results, list(tibble(
          approach = "mixture", distribution = dist, prefilter = pf,
          id = ss$id, condition = ss$condition,
          mean_rt = ss$mean_rt, var_rt = ss$var_rt
        )))
      }
    }
  }

  # --- Compute bias relative to clean ground truth ---
  all_results <- bind_rows(results) |>
    left_join(
      sumstats_clean |>
        select(id, condition, mean_rt, var_rt) |>
        rename(mean_rt_clean = mean_rt, var_rt_clean = var_rt),
      by = c("id", "condition")
    ) |>
    mutate(
      bias_mean_rt = mean_rt - mean_rt_clean,
      bias_var_rt  = var_rt - var_rt_clean
    )

  # Summarize per method-condition combination
  summary <- all_results |>
    summarize(
      mae_mean_rt = median(abs(bias_mean_rt), na.rm = TRUE),
      mae_var_rt  = median(abs(bias_var_rt), na.rm = TRUE),
      .by = c(approach, distribution, prefilter, condition)
    )

  # --- validate_fast_guesses on exgaussian flagging (no pre-filter) ---
  # Flag contaminants using the recommended default
  vfg_results <- tryCatch({
    flagged <- data_contam |>
      group_by(id, condition) |>
      mutate(
        contam_prob = flag_contaminant_rts(rt, distribution = "exgaussian")
      ) |>
      ungroup()

    # Apply validate_fast_guesses per participant x condition
    vfg <- flagged |>
      group_by(id, condition) |>
      reframe(
        validate_fast_guesses(
          contam_flag = contam_prob > 0.5,
          rt_data = rt,
          response = response
        ) |> as_tibble()
      )

    # Aggregate across participants
    c(
      bf_01_median = median(vfg$bf_01, na.rm = TRUE),
      pct_guess_in_hdi = mean(vfg$guess_in_hdi, na.rm = TRUE) * 100,
      pct_strong_evidence = mean(
        vfg$bf_evidence == "strong_for_guessing", na.rm = TRUE
      ) * 100
    )
  }, error = function(e) {
    c(bf_01_median = NA_real_,
      pct_guess_in_hdi = NA_real_,
      pct_strong_evidence = NA_real_)
  })

  # --- Return as a named vector for SimDesign ---
  out <- summary |>
    mutate(
      key = paste(approach, distribution, prefilter, condition, sep = "_")
    )

  ret <- setNames(out$mae_var_rt, paste0("mae_var_rt_", out$key))
  ret2 <- setNames(out$mae_mean_rt, paste0("mae_mean_rt_", out$key))
  c(ret, ret2, vfg_results)
}

###############################################################################!
# 4) Summarise Function --------------------------------------------------------
###############################################################################!

Summarise <- function(condition, results, fixed_objects) {
  # results is a matrix: rows = replications, cols = named statistics
  # SimDesign computes bias, RMSE, MAE etc. automatically if requested,
  # but here we just want the median and SD across replications

  ret <- apply(results, 2, function(x) {
    c(median = median(x, na.rm = TRUE),
      mean   = mean(x, na.rm = TRUE),
      sd     = sd(x, na.rm = TRUE))
  })

  # Flatten into a named vector
  as.vector(ret) |>
    setNames(paste0(rep(c("median_", "mean_", "sd_"),
                        each = 1, times = ncol(results)),
                    rep(colnames(results), each = 3)))
}

###############################################################################!
# 5) Run Simulation ------------------------------------------------------------
###############################################################################!

n_reps <- 500
n_cores <- min(parallel::detectCores() - 1, 15)

cat("=== Simulation Design ===\n")
cat("Conditions:", nrow(Design), "\n")
cat("Replications per condition:", n_reps, "\n")
cat("Total datasets:", nrow(Design) * n_reps, "\n")
cat("Parallel cores:", n_cores, "\n\n")

sim_results <- runSimulation(
  design       = Design,
  replications = n_reps,
  generate     = Generate,
  analyse      = Analyse,
  summarise    = Summarise,
  fixed_objects = fo,
  packages     = c("bmm", "tidyverse"),
  parallel     = TRUE,
  ncores       = n_cores,
  save         = TRUE,
  save_results = TRUE,
  filename     = "output/simulation_contaminant_methods",
  seed         = rep(2025, nrow(Design)),
  progress     = TRUE
)

###############################################################################!
# 6) Analyze & Visualize Results -----------------------------------------------
###############################################################################!

# Reshape SimDesign results for both var_rt and mean_rt
mae_cols <- grep("^median_mae_(var_rt|mean_rt)_", names(sim_results),
                 value = TRUE)

plot_data <- sim_results |>
  select(p_fastguess, p_slowcont, all_of(mae_cols)) |>
  pivot_longer(
    cols = all_of(mae_cols),
    names_to = "key",
    values_to = "median_mae"
  ) |>
  mutate(
    statistic = if_else(str_detect(key, "mae_var_rt"), "Var RT", "Mean RT"),
    key = str_remove(key, "^median_mae_(var_rt|mean_rt)_"),
    parts = str_split(key, "_"),
    approach = map_chr(parts, 1),
    distribution = map_chr(parts, 2),
    prefilter = map_chr(parts, 3),
    condition = map_chr(parts, 4),
    p_total = p_fastguess + p_slowcont,
    method_label = case_when(
      approach == "simple" ~ "Simple",
      approach == "robust" ~ "Robust",
      approach == "mixture" ~ paste0("Mixture (", distribution, ")")
    )
  ) |>
  select(-key, -parts)

## 6.1) Figure 1: Method comparison — 2x2 panel (statistic x condition) --------

# No-prefilter variants, 5 methods
method_levels <- c("Simple", "Robust",
                   "Mixture (exgaussian)", "Mixture (lognormal)",
                   "Mixture (invgaussian)")
method_labels <- c("Simple", "Robust",
                   "Mixture (ex-Gaussian)", "Mixture (lognormal)",
                   "Mixture (inv-Gaussian)")
method_colors <- c(
  "Simple" = "grey50",
  "Robust" = eam_palette[2],
  "Mixture (ex-Gaussian)" = eam_palette[1],
  "Mixture (lognormal)" = eam_palette[3],
  "Mixture (inv-Gaussian)" = eam_palette[4]
)
method_widths <- c(
  "Simple" = 0.6, "Robust" = 0.8,
  "Mixture (ex-Gaussian)" = 1.2,
  "Mixture (lognormal)" = 0.8,
  "Mixture (inv-Gaussian)" = 0.8
)

fig1_data <- plot_data |>
  filter(prefilter == "none") |>
  summarize(
    median_mae = median(median_mae, na.rm = TRUE),
    .by = c(method_label, condition, p_total, statistic)
  ) |>
  mutate(
    method_label = factor(method_label,
                          levels = method_levels, labels = method_labels),
    statistic = factor(statistic, levels = c("Mean RT", "Var RT")),
    condition = factor(condition, levels = c("easy", "hard"),
                       labels = c("Easy", "Hard"))
  )

p_fig1 <- fig1_data |>
  ggplot(aes(x = p_total, y = median_mae,
             color = method_label, group = method_label,
             linewidth = method_label)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(statistic ~ condition, scales = "free_y") +
  scale_color_manual(values = method_colors) +
  scale_linewidth_manual(values = method_widths, guide = "none") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Total Contamination Rate",
       y = "Median MAE",
       color = "Method") +
  clean_plot() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 2))

ggsave(here("figures", "02_method_comparison.pdf"),
       p_fig1, width = 6.5, height = 5, device = "pdf")

## 6.2) validate_fast_guesses summary ------------------------------------------

# Extract BF results from simulation
vfg_cols <- grep("^median_(bf_01|pct_guess|pct_strong)", names(sim_results),
                 value = TRUE)
if (length(vfg_cols) > 0) {
  cat("\n=== validate_fast_guesses Results ===\n\n")
  sim_results |>
    select(p_fastguess, p_slowcont, all_of(vfg_cols)) |>
    mutate(across(where(is.numeric), \(x) round(x, 2))) |>
    print(n = Inf)
}

###############################################################################!
# 7) Summary -------------------------------------------------------------------
###############################################################################!

cat("\n=== SIMULATION SUMMARY ===\n\n")

cat("Overall method ranking (no pre-filter, median MAE var_rt):\n\n")
plot_data |>
  filter(prefilter == "none") |>
  summarize(
    overall_mae = median(median_mae, na.rm = TRUE),
    .by = method_label
  ) |>
  arrange(overall_mae) |>
  mutate(overall_mae = round(overall_mae, 4)) |>
  print(n = Inf)

cat("\nOvercorrection at p = 0 (no contaminants):\n\n")
plot_data |>
  filter(prefilter == "none", p_fastguess == 0, p_slowcont == 0) |>
  select(method_label, condition, median_mae) |>
  pivot_wider(names_from = condition, values_from = median_mae) |>
  arrange(easy) |>
  mutate(across(where(is.numeric), \(x) round(x, 5))) |>
  print(n = Inf)

###############################################################################!
# 8) ndt x Pre-filter Interaction Analysis -------------------------------------
###############################################################################!

# This analysis shows that hard pre-filtering hurts the mixture method,
# especially at short ndt where the 150ms cutoff is close to legitimate RTs.

cat("\n=== Running ndt x pre-filter analysis ===\n")

ndt_values <- c(0.15, 0.20, 0.25, 0.30, 0.40)
n_reps_ndt <- 100
n_subjects_ndt <- 30
n_trials_ndt <- 200

ndt_results <- expand_grid(
  ndt_pop = ndt_values,
  p_fast = c(0, 0.05),
  p_slow = c(0, 0.05)
) |>
  pmap(function(ndt_pop, p_fast, p_slow) {
    cat("  ndt =", ndt_pop, "p_fast =", p_fast, "p_slow =", p_slow, "\n")
    map(1:n_reps_ndt, function(rep) {
      subj <- tibble(id = 1:n_subjects_ndt) |>
        mutate(
          drift = 2.0 + rnorm(n(), 0, fo$tau_drift),
          bound = 1.5 * exp(rnorm(n(), 0, fo$tau_bound)),
          ndt = pmax(0.05, ndt_pop + rnorm(n(), 0, fo$tau_ndt))
        )

      data_clean <- subj |>
        rowwise() |>
        reframe(
          id = id, trial = 1:n_trials_ndt,
          rddm(n_trials_ndt, drift = drift, bound = bound,
                ndt = ndt, zr = 0.5),
          ndt_true = ndt
        )

      ss_clean <- data_clean |>
        summarize(
          ezdm_summary_stats(rt, response, version = "3par",
                             method = "simple"),
          .by = id
        )

      # Inject contaminants
      if (p_fast > 0 || p_slow > 0) {
        rt_stats <- data_clean |>
          group_by(id) |>
          summarize(m = mean(rt), s = sd(rt), .groups = "drop")
        data_contam <- data_clean |>
          left_join(rt_stats, by = "id") |>
          group_by(id) |>
          mutate(
            n = n(), draw = sample(n()),
            is_fast = draw <= round(n * p_fast),
            is_slow = !is_fast & draw > (n - round(n * p_slow))
          ) |>
          ungroup() |>
          mutate(
            rt = if_else(is_fast, runif(n(), 0.05, ndt_true), rt),
            response = if_else(is_fast,
                               as.numeric(rbinom(n(), 1, 0.5)), response),
            rt = if_else(is_slow,
                         runif(n(), m + 2*s, m + 5*s), rt)
          )
      } else {
        data_contam <- data_clean
      }

      data_hard <- data_contam |>
        filter(rt >= fo$hard_lower & rt <= fo$hard_upper)

      eval_method <- function(d, method, dist = "exgaussian") {
        tryCatch({
          ss <- d |>
            summarize(
              ezdm_summary_stats(rt, response, version = "3par",
                                 method = method, distribution = dist),
              .by = id
            ) |>
            left_join(ss_clean |> select(id, mean_rt, var_rt) |>
                        rename(mean_rt_clean = mean_rt,
                               var_rt_clean = var_rt),
                      by = "id")
          tibble(
            mae_mean_rt = median(abs(ss$mean_rt - ss$mean_rt_clean),
                                  na.rm = TRUE),
            mae_var_rt = median(abs(ss$var_rt - ss$var_rt_clean),
                                 na.rm = TRUE)
          )
        }, error = function(e) tibble(mae_mean_rt = NA, mae_var_rt = NA))
      }

      bind_rows(
        eval_method(data_contam, "simple") |>
          mutate(method = "Simple", prefilter = "No filter"),
        eval_method(data_hard, "simple") |>
          mutate(method = "Simple", prefilter = "Hard filter"),
        eval_method(data_contam, "robust") |>
          mutate(method = "Robust", prefilter = "No filter"),
        eval_method(data_hard, "robust") |>
          mutate(method = "Robust", prefilter = "Hard filter"),
        eval_method(data_contam, "mixture") |>
          mutate(method = "Mixture (ex-Gaussian)", prefilter = "No filter"),
        eval_method(data_hard, "mixture") |>
          mutate(method = "Mixture (ex-Gaussian)", prefilter = "Hard filter")
      ) |>
        mutate(ndt_pop = ndt_pop, p_fast = p_fast, p_slow = p_slow,
               rep = rep)
    }) |> bind_rows()
  }) |> bind_rows()

# Save ndt results
write_csv(ndt_results, here("data", "simulation_ndt_prefilter.csv"))

## 8.1) Figure 2: Pre-filter effect by ndt (two-panel: mean_rt + var_rt) -------

fig2_data <- ndt_results |>
  filter(p_fast == 0.05, p_slow == 0.05) |>
  summarize(
    mae_mean_rt = median(mae_mean_rt, na.rm = TRUE),
    mae_var_rt = median(mae_var_rt, na.rm = TRUE),
    .by = c(ndt_pop, method, prefilter)
  ) |>
  pivot_longer(
    cols = c(mae_mean_rt, mae_var_rt),
    names_to = "statistic",
    values_to = "mae"
  ) |>
  mutate(
    statistic = factor(statistic,
                       levels = c("mae_mean_rt", "mae_var_rt"),
                       labels = c("Mean RT", "Var RT"))
  )

ndt_colors <- c(
  "Simple" = "grey50",
  "Robust" = eam_palette[2],
  "Mixture (ex-Gaussian)" = eam_palette[1]
)

p_fig2 <- fig2_data |>
  ggplot(aes(x = ndt_pop, y = mae,
             color = method, linetype = prefilter,
             group = interaction(method, prefilter))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ statistic, scales = "free_y") +
  scale_color_manual(values = ndt_colors) +
  scale_linetype_manual(values = c(
    "No filter" = "solid",
    "Hard filter" = "dashed"
  )) +
  scale_x_continuous(breaks = ndt_values,
                     labels = paste0(ndt_values * 1000, " ms")) +
  labs(x = "Population Non-Decision Time",
       y = "Median MAE",
       color = "Method",
       linetype = "Pre-filter") +
  clean_plot() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8))

ggsave(here("figures", "02_prefilter_ndt.pdf"),
       p_fig2, width = 6.5, height = 3.5, device = "pdf")

cat("\n=== ndt x pre-filter analysis complete ===\n")
ndt_results |>
  filter(p_fast == 0.05, p_slow == 0.05) |>
  summarize(
    mae_mean_rt = median(mae_mean_rt, na.rm = TRUE),
    mae_var_rt = median(mae_var_rt, na.rm = TRUE),
    .by = c(ndt_pop, method, prefilter)
  ) |>
  pivot_wider(names_from = prefilter,
              values_from = c(mae_mean_rt, mae_var_rt)) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  print(n = Inf)
