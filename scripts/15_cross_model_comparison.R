#' Cross-Model Comparison
#'
#' Loads all fitted models from Tutorials 1-4 and compares parameter estimates
#' across EZDM, CSWALD, and DDM for both accuracy-coded and response-coded
#' analyses. Overlays reference values from Wagenmakers et al. (2008), Table 2
#' (converted to s = 1 scaling).

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# -- Load all fitted models ----------------------------------------------------
fit_acc_ezdm      <- readRDS(here("output", "fit_ezdm_lexical_acc.rds"))
fit_acc_ezdm_freq <- readRDS(here("output", "fit_ezdm_lexical_acc_freq.rds"))
fit_cswald_simple <- readRDS(here("output", "fit_cswald_lexical_simple.rds"))
fit_ddm_acc       <- readRDS(here("output", "fit_ddm_lexical_acc.rds"))

fit_resp_ezdm     <- readRDS(here("output", "fit_ezdm_lexical_resp.rds"))
fit_resp_ezdm_freq <- readRDS(here("output", "fit_ezdm_lexical_resp_freq.rds"))
fit_cswald_crisk  <- readRDS(here("output", "fit_cswald_lexical_resp.rds"))
fit_ddm_resp      <- readRDS(here("output", "fit_ddm_lexical_resp.rds"))

# -- Load reference parameters (Wagenmakers et al., 2008, s = 1) ---------------
ref_params <- read_csv(
  here("data", "wagenmakers2008_reference_params.csv"),
  show_col_types = FALSE
)

###############################################################################!
# Helper: Extract fixed effects ------------------------------------------------
###############################################################################!

extract_fixef <- function(fit, model_name, params, transforms = list()) {
  fe <- fixef(fit)
  map_dfr(params, function(par_name) {
    rows <- grep(paste0("^", par_name, "_"), rownames(fe), value = TRUE)
    map_dfr(rows, function(r) {
      est <- fe[r, ]
      tfm <- transforms[[par_name]] %||% identity
      tibble(
        model = model_name,
        parameter = par_name,
        level = sub(paste0("^", par_name, "_"), "", r),
        estimate = tfm(est["Estimate"]),
        lower = tfm(est["Q2.5"]),
        upper = tfm(est["Q97.5"])
      )
    })
  })
}

clean_levels <- function(df) {
  df |>
    mutate(
      level = str_remove(level, "^(stim_cat|condition|freq)"),
      level = str_remove(level, "Intercept")
    )
}

###############################################################################!
# 1) Accuracy-Coded Comparison -------------------------------------------------
###############################################################################!

# EZDM (stim_cat), EZDM freq, CSWALD (stim_cat), DDM (freq)
# Compare drift at both granularities, bound across all models

acc_transforms <- list(bound = exp)

# Bound comparison: all 4 models (CSWALD simple = half boundary)
bound_acc <- bind_rows(
  extract_fixef(fit_acc_ezdm, "EZDM", "bound", acc_transforms),
  extract_fixef(fit_acc_ezdm_freq, "EZDM (freq)", "bound", acc_transforms),
  extract_fixef(fit_cswald_simple, "CSWALD simple", "bound", acc_transforms),
  extract_fixef(fit_ddm_acc, "DDM", "bound", acc_transforms)
) |>
  clean_levels() |>
  mutate(
    model = factor(model, levels = c("EZDM", "EZDM (freq)", "CSWALD simple", "DDM")),
    parameter = "bound"
  )

# Reference values for bound
ref_bound <- ref_params |>
  filter(parameter == "bound") |>
  transmute(level = condition, ref_value = value_s1)

# Drift by stim_cat: EZDM and CSWALD
drift_stim <- bind_rows(
  extract_fixef(fit_acc_ezdm, "EZDM", "drift"),
  extract_fixef(fit_cswald_simple, "CSWALD simple", "drift")
) |> clean_levels()

# Drift by freq: EZDM freq and DDM
drift_freq <- bind_rows(
  extract_fixef(fit_acc_ezdm_freq, "EZDM (freq)", "drift"),
  extract_fixef(fit_ddm_acc, "DDM", "drift")
) |> clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM (freq)", "DDM")))

# Reference values for drift (accuracy-coded: all positive, NW = |NW|)
ref_drift_freq <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(
    level = level,
    ref_value = ifelse(level == "nonword", abs(value_s1), value_s1)
  )

# Plot drift by frequency: EZDM freq vs DDM with reference
p_drift_freq <- ggplot(drift_freq, aes(x = level, y = estimate,
                                        color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  geom_point(data = ref_drift_freq, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 4)]) +
  labs(x = NULL, y = "Estimate",
       title = "Drift by frequency",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot drift by stim_cat: EZDM vs CSWALD
p_drift_stim <- ggplot(drift_stim, aes(x = level, y = estimate,
                                         color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.4), size = 0.5) +
  scale_color_manual(values = eam_palette[c(1, 3)]) +
  labs(x = NULL, y = "Estimate",
       title = "Drift by stim_cat",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot bound comparison with reference
p_bound <- ggplot(bound_acc, aes(x = level, y = estimate,
                                  color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(data = ref_bound, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 2, 3, 4)]) +
  labs(x = NULL, y = "Estimate",
       title = "Bound by condition",
       color = "Model", shape = "Model") +
  annotate("text", x = 0.5, y = Inf, label = "x = Wagenmakers (2008)",
           hjust = 0, vjust = 1.5, size = 2.5) +
  clean_plot()

p_acc <- (p_drift_freq | p_drift_stim) / p_bound +
  plot_layout(heights = c(1, 1))

ggsave(here("figures", "15_cross_model_accuracy.pdf"),
       p_acc, width = 6.5, height = 6)

###############################################################################!
# 2) Response-Coded Comparison -------------------------------------------------
###############################################################################!

resp_transforms <- list(bound = exp, zr = plogis)

# Bound: all 4 models
bound_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm, "EZDM 4par", "bound", resp_transforms),
  extract_fixef(fit_resp_ezdm_freq, "EZDM 4par (freq)", "bound", resp_transforms),
  extract_fixef(fit_cswald_crisk, "CSWALD crisk", "bound", resp_transforms),
  extract_fixef(fit_ddm_resp, "DDM", "bound", resp_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM 4par", "EZDM 4par (freq)",
                                           "CSWALD crisk", "DDM")))

# zr: all 4 models
zr_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm, "EZDM 4par", "zr", resp_transforms),
  extract_fixef(fit_resp_ezdm_freq, "EZDM 4par (freq)", "zr", resp_transforms),
  extract_fixef(fit_cswald_crisk, "CSWALD crisk", "zr", resp_transforms),
  extract_fixef(fit_ddm_resp, "DDM", "zr", resp_transforms)
) |>
  clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM 4par", "EZDM 4par (freq)",
                                           "CSWALD crisk", "DDM")))

# Reference values for zr
ref_zr <- ref_params |>
  filter(parameter == "zr") |>
  transmute(level = condition, ref_value = value_s1)

# Drift by stim_cat: EZDM 4par and CSWALD crisk
drift_stim_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm, "EZDM 4par", "drift"),
  extract_fixef(fit_cswald_crisk, "CSWALD crisk", "drift")
) |> clean_levels()

# Drift by freq: EZDM 4par freq and DDM
drift_freq_resp <- bind_rows(
  extract_fixef(fit_resp_ezdm_freq, "EZDM 4par (freq)", "drift"),
  extract_fixef(fit_ddm_resp, "DDM", "drift")
) |> clean_levels() |>
  mutate(model = factor(model, levels = c("EZDM 4par (freq)", "DDM")))

ref_drift_resp <- ref_params |>
  filter(parameter == "drift", level != "word_avg") |>
  transmute(level = level, ref_value = value_s1)

# Plot response-coded drift by frequency: EZDM freq vs DDM
p_drift_resp <- ggplot(drift_freq_resp, aes(x = level, y = estimate,
                                             color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.5), size = 0.5) +
  geom_point(data = ref_drift_resp, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = eam_palette[c(1, 4)]) +
  labs(x = NULL, y = "Estimate", title = "Drift by frequency",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot response-coded drift by stim_cat: EZDM vs CSWALD
p_drift_stim_resp <- ggplot(
  drift_stim_resp,
  aes(x = level, y = estimate, color = model, shape = model)
) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.4), size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = eam_palette[c(1, 3)]) +
  labs(x = NULL, y = "Estimate",
       title = "Drift by stim_cat",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot bound with reference
p_bound_resp <- ggplot(bound_resp, aes(x = level, y = estimate,
                                        color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(data = ref_bound, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  scale_color_manual(values = eam_palette[c(1, 2, 3, 4)]) +
  labs(x = NULL, y = "Estimate", title = "Bound by condition",
       color = "Model", shape = "Model") +
  clean_plot()

# Plot zr with reference
p_zr <- ggplot(zr_resp, aes(x = level, y = estimate,
                              color = model, shape = model)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(data = ref_zr, aes(x = level, y = ref_value),
             inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = eam_palette[c(1, 2, 3, 4)]) +
  labs(x = NULL, y = "Estimate", title = "Rel. starting point by condition",
       color = "Model", shape = "Model") +
  clean_plot()

p_resp <- (p_drift_resp | p_drift_stim_resp) /
  (p_bound_resp | p_zr) +
  plot_layout(heights = c(1, 1))

ggsave(here("figures", "15_cross_model_response.pdf"),
       p_resp, width = 6.5, height = 6)

###############################################################################!
# 3) Computational Cost Summary ------------------------------------------------
###############################################################################!

cat("\n=== Computational Cost Summary ===\n")
cat("These times are approximate and depend on hardware:\n")
cat("  EZDM 3par (accuracy):    ~10 seconds\n")
cat("  EZDM 4par (response):    ~15 seconds\n")
cat("  CSWALD simple (accuracy): ~15 minutes\n")
cat("  CSWALD crisk (response):  ~30 minutes\n")
cat("  DDM accuracy:             ~60 minutes\n")
cat("  DDM response:             ~90 minutes\n")
