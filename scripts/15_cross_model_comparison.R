# Part III: Cross-approach comparison at the freq specification
# Fits the EZ-Diffusion Model with freq on drift, checks the trial-level fits,
# and compares the EZ-Diffusion, Censored-Shifted Wald, and DDM estimates for
# both codings against Wagenmakers et al. (2008, Table 2).

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
# Every fit uses 2 threads per chain so that fitting times are comparable
chains <- 5
cores  <- 5
threads <- 2
warmup <- 1000
iter   <- 2000
adapt_delta   <- 0.90
max_treedepth <- 12

# -- Load preprocessed data (created by scripts/01_prepare_lexical_data.R) ----
freq_levels <- c("high", "low", "very_low", "nonword")

data_accuracy_clean <- read_csv(
  here("data", "lexical_accuracy_triallevel.csv"), show_col_types = FALSE
) |>
  mutate(freq = factor(freq, levels = freq_levels))

data_response_clean <- read_csv(
  here("data", "lexical_response_triallevel.csv"), show_col_types = FALSE
) |>
  mutate(freq = factor(freq, levels = freq_levels))

# -- Load fitted models from scripts 11-14 -------------------------------------
fit_acc_ezdm      <- readRDS(here("output", "fit_ezdm_lexical_acc.rds"))
fit_cswald_simple <- readRDS(here("output", "fit_cswald_lexical_simple.rds"))
fit_ddm_acc       <- readRDS(here("output", "fit_ddm_lexical_acc.rds"))

fit_resp_ezdm     <- readRDS(here("output", "fit_ezdm_lexical_resp.rds"))
fit_cswald_crisk  <- readRDS(here("output", "fit_cswald_lexical_resp.rds"))
fit_ddm_resp      <- readRDS(here("output", "fit_ddm_lexical_resp.rds"))

# -- Reference parameters (Wagenmakers et al., 2008, rescaled to s = 1) -------
ref_params <- read_csv(
  here("data", "wagenmakers2008_reference_params.csv"),
  show_col_types = FALSE
)

###############################################################################!
# 1) EZ-Diffusion Fits at the freq Specification ------------------------------
###############################################################################!

## 1.1) EZ-Diffusion 3par with freq (accuracy-coded) ---------------------------

# The frequency refit of Tutorial 1: summary statistics per id x condition x freq
data_acc_ezdm_freq <- data_accuracy_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "3par",
                       method = "mixture", distribution = "exgaussian"),
    .by = c(id, condition, freq)
  )

model_ezdm <- ezdm(
  mean_rt  = "mean_rt",
  var_rt   = "var_rt",
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "3par"
)

formula_ezdm_freq <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

fit_acc_ezdm_freq <- bmm(
  formula = formula_ezdm_freq,
  data    = data_acc_ezdm_freq,
  model   = model_ezdm,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_acc_freq")
)

## 1.2) EZ-Diffusion 4par with freq (response-coded) --------------------------

# Robust moments with min_trials = 3 keep cells with few minority-boundary responses
data_resp_ezdm_freq <- data_response_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par",
                       method = "robust",
                       min_trials = 3),
    .by = c(id, condition, freq)
  )

# Trial counts on the minority boundary
minority_counts <- data_resp_ezdm_freq |>
  mutate(n_lower = n_trials - n_upper,
         n_minority = pmin(n_upper, n_lower)) |>
  group_by(condition, freq) |>
  summarise(
    median_minority = median(n_minority),
    min_minority = min(n_minority),
    .groups = "drop"
  )
cat("\n=== EZ-Diffusion 4par: Trial counts on minority boundary ===\n")
print(minority_counts)

model_ezdm_4par <- ezdm(
  mean_rt  = c("mean_rt_upper", "mean_rt_lower"),
  var_rt   = c("var_rt_upper", "var_rt_lower"),
  n_upper  = "n_upper",
  n_trials = "n_trials",
  version  = "4par"
)

formula_ezdm_resp_freq <- bmf(
  drift ~ 0 + freq + (0 + freq | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + condition + (0 + condition | id)
)

fit_resp_ezdm_freq <- bmm(
  formula = formula_ezdm_resp_freq,
  data    = data_resp_ezdm_freq,
  model   = model_ezdm_4par,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_resp_freq")
)

## 1.3) The same 4par fit from mixture moments ---------------------------------

# Isolates the offset of the robust moments: same cells and formula, mixture
# moments wherever the EM converges (robust fallback otherwise, counted below)
data_resp_ezdm_freq_mix <- data_response_clean |>
  summarize(
    ezdm_summary_stats(rt, response, version = "4par",
                       method = "mixture", distribution = "exgaussian",
                       min_trials = 3),
    .by = c(id, condition, freq)
  )
mix_fallback <- data_resp_ezdm_freq_mix |>
  summarise(
    cells = n(),
    boundaries_fallback = sum(is.na(contaminant_prop_upper) & !is.na(mean_rt_upper)) +
                          sum(is.na(contaminant_prop_lower) & !is.na(mean_rt_lower)),
    cells_fallback = sum((is.na(contaminant_prop_upper) & !is.na(mean_rt_upper)) |
                         (is.na(contaminant_prop_lower) & !is.na(mean_rt_lower)))
  )
cat("\n=== EZ-Diffusion 4par, mixture moments: EM fallback to robust ===\n")
print(mix_fallback)
write_csv(mix_fallback, here("data", "lexical_response_ezdm4par_freq_mix_fallback.csv"))

fit_resp_ezdm_freq_mix <- bmm(
  formula = formula_ezdm_resp_freq,
  data    = data_resp_ezdm_freq_mix,
  model   = model_ezdm_4par,
  chains  = chains,
  cores   = cores,
  threads = threading(threads),
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_ezdm_lexical_resp_freq_mix")
)

###############################################################################!
# 2) Composite Posterior Predictive Check (Part II figure) ---------------------
###############################################################################!

# One row per trial-level fit; the signed RT (lower-boundary responses
# negative) checks choices and RTs jointly. Each row uses one set of 50 draws.
trial_fits <- list(
  "Censored-Shifted Wald (simple), accuracy-coded" = fit_cswald_simple,
  "Full DDM, accuracy-coded"                       = fit_ddm_acc,
  "Censored-Shifted Wald (crisk), response-coded"  = fit_cswald_crisk,
  "Full DDM, response-coded"                       = fit_ddm_resp
)
trial_vars <- c(rt = "Response time", response = "Response (0 = lower, 1 = upper)",
                signed_rt = "Signed response time")
set.seed(1)
p_pp_trial <- imap(trial_fits, \(fit, model_lbl) {
  draws <- sample(ndraws(fit), 50)
  panels <- imap(trial_vars, \(var_lbl, v)
    pp_check(fit, resp_var = v, draw_ids = draws) +
      labs(title = var_lbl, subtitle = NULL) + clean_plot() +
      theme(legend.position = "none", plot.title = element_text(size = 8)) +
      # Observed RTs lie below 3 s; long simulated tails would flatten the density
      coord_cartesian(xlim = switch(v, rt = c(0, 3), signed_rt = c(-3, 3), NULL))
  )
  wrap_elements(full = grid::textGrob(model_lbl, x = 0, hjust = 0,
                                      gp = grid::gpar(fontsize = 9, fontface = "bold"))) /
    wrap_plots(panels, nrow = 1) +
    plot_layout(heights = c(0.08, 1))
}) |>
  wrap_plots(ncol = 1)
ggsave(here("figures", "15_pp_check_trial_level.pdf"),
       p_pp_trial, width = 6.5, height = 5)

## 2.1) Per-cell interval checks for the trial-level fits ---------------------

# Counts the id x condition x freq cells whose mean RT, RT variance, or proportion
# upper fall outside the 90% predictive interval, as done for the EZ-Diffusion fits
cell_check <- function(fit, draws) {
  d <- fit$data
  resp_col <- intersect(c("dec", "response"), names(d))[1]
  resp <- d[[resp_col]]
  resp <- if (is.numeric(resp)) resp else as.numeric(resp == "upper")
  cell <- interaction(d$id, d$condition, d$freq, drop = TRUE)
  cell_stats <- function(rt, upper) {
    tibble(cell = cell, rt = rt, upper = upper) |>
      summarise(mean_rt = mean(rt), var_rt = var(rt), p_upper = mean(upper),
                .by = cell)
  }
  obs <- cell_stats(d$rt, resp)
  # pp_simulate() returns ndraws x nobs matrices of rt and response; it is much
  # faster than posterior_predict() for these models
  prep <- brms::prepare_predictions(fit, draw_ids = draws)
  yrep <- bmm::pp_simulate(fit$bmm$model, prep)
  sims <- map(seq_len(nrow(yrep$rt)),
              \(i) cell_stats(yrep$rt[i, ], as.numeric(yrep$response[i, ]))) |>
    list_rbind(names_to = "draw")
  lims <- sims |>
    summarise(across(c(mean_rt, var_rt, p_upper),
                     list(lo = \(x) quantile(x, 0.05), hi = \(x) quantile(x, 0.95))),
              .by = cell)
  obs |>
    left_join(lims, join_by(cell)) |>
    summarise(cells   = n(),
              mean_rt = sum(mean_rt < mean_rt_lo | mean_rt > mean_rt_hi),
              var_rt  = sum(var_rt  < var_rt_lo  | var_rt  > var_rt_hi),
              p_upper = sum(p_upper < p_upper_lo | p_upper > p_upper_hi))
}
set.seed(1)
trial_cell_checks <- imap(trial_fits, \(fit, lbl) {
  draws <- sample(ndraws(fit), 50)
  cell_check(fit, draws) |> mutate(fit = lbl, .before = 1)
}) |>
  list_rbind()
cat("\n=== Trial-level fits: cells outside their 90% predictive intervals ===\n")
print(trial_cell_checks)
write_csv(trial_cell_checks, here("data", "trial_level_ppc_cell_counts.csv"))

###############################################################################!
# 3) Accuracy-Coded Comparison -------------------------------------------------
###############################################################################!

# Drift by frequency, bound by condition, and ndt; link functions are inverted
acc_transforms <- list(bound = exp, ndt = exp)

# The simple version's bound is half the boundary separation (so x 2), and its
# drift has a log link
cswald_transforms <- list(drift = exp, bound = \(x) exp(x) * 2, ndt = exp)

acc_levels <- c("EZ-Diffusion", "Censored-Shifted Wald (simple)", "Full DDM")
acc_all <- bind_rows(
  extract_fixef(fit_acc_ezdm_freq, acc_levels[1], c("drift", "bound", "ndt"), acc_transforms),
  extract_fixef(fit_cswald_simple, acc_levels[2], c("drift", "bound", "ndt"), cswald_transforms),
  extract_fixef(fit_ddm_acc,       acc_levels[3], c("drift", "bound", "ndt"), acc_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = acc_levels),
         level = if_else(level == "", "ndt", level))

# Reference values; with accuracy coding the nonword drift is positive
ref_drift_freq <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(level = level,
            ref_value = ifelse(level == "nonword", abs(value_s1), value_s1))
ref_bound <- ref_params |>
  filter(parameter == "bound") |>
  transmute(level = condition, ref_value = value_s1)
ref_ndt <- ref_params |>
  filter(parameter == "ndt") |>
  transmute(level = "ndt", ref_value = value_s1)

cmp_panel <- function(d, ref, title, palette, y_lab = NULL, dodge = 0.6) {
  ggplot(d, aes(x = level, y = estimate, color = model, shape = model)) +
    geom_pointrange(aes(ymin = lower, ymax = upper),
                    position = position_dodge(width = dodge), size = 0.5) +
    geom_point(data = ref, aes(x = level, y = ref_value),
               inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
    scale_color_manual(values = palette) +
    labs(x = NULL, y = y_lab, title = title, color = "Model", shape = "Model") +
    clean_plot() +
    theme(plot.title = element_text(size = 9))
}

p_acc_drift <- cmp_panel(filter(acc_all, parameter == "drift") |>
                           mutate(level = factor(level, levels = freq_levels)),
                         ref_drift_freq, "Drift",
                         eam_palette[c(1, 3, 4)], y_lab = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_acc_bound <- cmp_panel(filter(acc_all, parameter == "bound"), ref_bound,
                         "Bound", eam_palette[c(1, 3, 4)])
p_acc_ndt   <- cmp_panel(filter(acc_all, parameter == "ndt"), ref_ndt,
                         "ndt", eam_palette[c(1, 3, 4)])

p_acc <- (p_acc_drift | p_acc_bound | p_acc_ndt) +
  plot_layout(guides = "collect", widths = c(2, 1.2, 0.8)) &
  theme(legend.position = "bottom")

ggsave(here("figures", "15_cross_model_accuracy.pdf"),
       p_acc, width = 6.5, height = 3.8)

###############################################################################!
# 4) Response-Coded Comparison -------------------------------------------------
###############################################################################!

# Drift by frequency, bound by condition, and zr by condition
resp_transforms <- list(bound = exp, zr = plogis)
resp_levels <- c("EZ-Diffusion (4par)", "Censored-Shifted Wald (crisk)", "Full DDM")
resp_all <- bind_rows(
  extract_fixef(fit_resp_ezdm_freq, resp_levels[1], c("drift", "bound", "zr"), resp_transforms),
  extract_fixef(fit_cswald_crisk,   resp_levels[2], c("drift", "bound", "zr"), resp_transforms),
  extract_fixef(fit_ddm_resp,       resp_levels[3], c("drift", "bound", "zr"), resp_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = resp_levels))

ref_drift_resp <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(level = level, ref_value = value_s1)
ref_zr <- ref_params |>
  filter(parameter == "zr") |>
  transmute(level = condition, ref_value = value_s1)

p_resp_drift <- cmp_panel(filter(resp_all, parameter == "drift") |>
                            mutate(level = factor(level, levels = freq_levels)),
                          ref_drift_resp, "Drift",
                          eam_palette[c(1, 3, 4)], y_lab = "Estimate") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_resp_bound <- cmp_panel(filter(resp_all, parameter == "bound"), ref_bound,
                          "Bound", eam_palette[c(1, 3, 4)])
p_resp_zr <- cmp_panel(filter(resp_all, parameter == "zr"), ref_zr,
                       "Starting point", eam_palette[c(1, 3, 4)]) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.3)

p_resp <- (p_resp_drift | p_resp_bound | p_resp_zr) +
  plot_layout(guides = "collect", widths = c(2, 1.2, 1)) &
  theme(legend.position = "bottom")

ggsave(here("figures", "15_cross_model_response.pdf"),
       p_resp, width = 6.5, height = 3.8)

###############################################################################!
# 5) Computational Cost Summary ------------------------------------------------
###############################################################################!

# Wall-clock time of the slowest chain, read from the cached fits
cost <- tibble(
  fit = c("EZ-Diffusion 3par (stim_cat)", "EZ-Diffusion 3par (freq)", "EZ-Diffusion 4par (stim_cat)",
          "EZ-Diffusion 4par (freq)", "Censored-Shifted Wald simple", "Censored-Shifted Wald crisk",
          "DDM accuracy", "DDM response"),
  hours = purrr::map_dbl(
    list(fit_acc_ezdm, fit_acc_ezdm_freq, fit_resp_ezdm, fit_resp_ezdm_freq,
         fit_cswald_simple, fit_cswald_crisk, fit_ddm_acc, fit_ddm_resp),
    \(f) max(rowSums(rstan::get_elapsed_time(f$fit))) / 3600
  )
)
cat("\n=== Computational Cost Summary (wall-clock, slowest chain) ===\n")
print(cost |> mutate(hours = round(hours, 3)), n = Inf)
