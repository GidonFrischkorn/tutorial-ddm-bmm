#' Tutorial 2: EZ-Diffusion Model (4-Parameter) with bmm
#'
#' This tutorial extends the EZDM to four parameters by adding the relative
#' starting point (zr). The 4par version requires separate summary statistics
#' for upper and lower boundary responses. We also compare the 3par and 4par
#' models using bridge sampling.

###############################################################################!
# 0) R Setup: Packages & Data -------------------------------------------------
###############################################################################!

pacman::p_load(
  here, bmm, brms, bridgesampling,
  tidyverse, tidybayes, patchwork
)
source(here("functions", "clean_plot.R"))

# -- Fitting settings ----------------------------------------------------------
chains <- 4
cores  <- 4
warmup <- 1000
iter   <- 2000

# -- Load aggregated data (4par version) ---------------------------------------
data_ezdm_4par <- read_csv(here("data", "simulated_ezdm_4par.csv"))

###############################################################################!
# 1) Model Specification -------------------------------------------------------
###############################################################################!

# TODO: Specify the 4par EZDM model
#
# model <- ezdm(
#   resp_vars = c("mean_rt", "var_rt", "n_upper", "n_trials"),
#   version = "4par"
# )
#
# formula <- bmf(
#   drift ~ 0 + condition + (0 + condition | id),
#   bound ~ 1 + (1 | id),
#   ndt   ~ 1 + (1 | id),
#   zr    ~ 1 + (1 | id)
# )

###############################################################################!
# 2) Model Fitting -------------------------------------------------------------
###############################################################################!

# TODO: Fit the 4par model
# TODO: Load 3par model from Tutorial 1 for comparison

###############################################################################!
# 3) Model Evaluation -----------------------------------------------------------
###############################################################################!

# TODO: Parameter estimates (focus on zr interpretation)
# TODO: Model comparison: 3par vs 4par (bridge sampling / LOO)
# TODO: When is the 4par version preferred?
