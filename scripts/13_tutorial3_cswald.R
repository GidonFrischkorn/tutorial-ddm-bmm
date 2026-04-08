#' Tutorial 3: Trial-Level Modeling with the CSWALD
#'
#' Fits both CSWALD versions to Wagenmakers et al. (2008) lexical decision data:
#' - Simple (censored): accuracy-coded, appropriate for low error rates (<20%)
#' - Competing-risks (crisk): response-coded, estimates zr, handles higher errors
#'
#' Demonstrates when to use each version and the known ~10% bound
#' underestimation in the crisk version.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms, rtdists,
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

# -- Load preprocessed data ---------------------------------------------------
# Run scripts/prepare_lexical_data.R once first to generate these files
data_accuracy_clean <- read_csv(here("data", "lexical_accuracy_triallevel.csv"),
                                 show_col_types = FALSE)
data_response_clean <- read_csv(here("data", "lexical_response_triallevel.csv"),
                                 show_col_types = FALSE)

###############################################################################!
# 1) CSWALD Simple: Accuracy-Coded Analysis ------------------------------------
###############################################################################!

model_simple <- cswald(
  rt       = "rt",
  response = "response",
  version  = "simple"
)

formula_simple <- bmf(
  drift ~ 0 + stim_cat + (0 + stim_cat | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id)
)

fit_cswald_simple <- bmm(
  formula = formula_simple,
  data    = data_accuracy_clean,
  model   = model_simple,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_cswald_lexical_simple")
)

## Convergence & estimates -----------------------------------------------------
summary(fit_cswald_simple)
rhat_simple <- max(brms::rhat(fit_cswald_simple), na.rm = TRUE)

## Posterior predictive check ---------------------------------------------------
p_pp_simple <- pp_check(fit_cswald_simple, ndraws = 50) + clean_plot()
ggsave(here("figures", "13_pp_check_simple.pdf"), p_pp_simple,
       width = 6.5, height = 4, device = "pdf")


###############################################################################!
# 2) CSWALD Competing-Risks: Response-Coded Analysis ---------------------------
###############################################################################!

model_crisk <- cswald(
  rt       = "rt",
  response = "response",
  version  = "crisk"
)

formula_crisk <- bmf(
  drift ~ 0 + stim_cat + (0 + stim_cat | id),
  bound ~ 0 + condition + (0 + condition | id),
  ndt   ~ 1 + (1 | id),
  zr    ~ 0 + condition + (0 + condition | id)
)

fit_cswald_crisk <- bmm(
  formula = formula_crisk,
  data    = data_response_clean,
  model   = model_crisk,
  chains  = chains,
  cores   = cores,
  warmup  = warmup,
  iter    = iter,
  backend = "cmdstanr",
  sample_prior = "yes",
  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
  file    = here("output", "fit_cswald_lexical_resp")
)

## Convergence & estimates -----------------------------------------------------
summary(fit_cswald_crisk)
rhat_crisk <- max(brms::rhat(fit_cswald_crisk), na.rm = TRUE)

## Posterior predictive check ---------------------------------------------------
p_pp_crisk <- pp_check(fit_cswald_crisk, ndraws = 50) + clean_plot()
ggsave(here("figures", "13_pp_check_crisk.pdf"), p_pp_crisk,
       width = 6.5, height = 4, device = "pdf")
