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
#' Design (v7): the mixed process is the primary result, the pure processes
#'         are mechanism diagnostics: 3 pure processes x 2 rates (.02, .10)
#'         + mixed at Beta mean rates .02/.05/.10 + shared clean baseline,
#'         x 2 trial counts = 20 rows x 200 replications, 50 subjects each.
#' Contaminant processes (contaminants replace round(n * p) trials per cell;
#' the pre-injection clean trials are returned alongside as the eCDF
#' calibration target — replacement means they are otherwise lost):
#'   leading_edge  anticipations straddling the leading edge: each subject
#'               draws a depth x_i ~ U(0.1, 0.5), shared across conditions;
#'               RT ~ U((1 - x_i) ndt_i, (1 + x_i) ndt_i), response at chance
#'               (Lerche, Voss & Nagler, 2017, adapted to centre on ndt_i) —
#'               straddles both the hard cutoffs and the core's leading edge,
#'               so neither the cutoffs nor the mixture get a free case
#'   delay       delayed start-up: U(0, 2 s) added to a genuine RT, response
#'               kept (Ratcliff & Tuerlinckx, 2002; Ratcliff, 2008;
#'               Wagenmakers et al., 2008) — partly overlaps the core
#'   guess       informationless response: the trial is redrawn from the DDM
#'               with drift = 0 (Vandekerckhove, Tuerlinckx & Lee, 2008;
#'               Chavez De La Pena, Shin & Vandekerckhove, 2026) — fully
#'               overlaps the core, accuracy at chance
#'   mixed       realistic cocktail: each subject draws a total rate
#'               p_i ~ Beta(mean p_contam, concentration 10; rates spread
#'               ~0-25% at mean .05, i.e. many near-clean subjects plus a
#'               minority of heavy contaminators) and a process mix
#'               w_i ~ Dirichlet(1, 1, 1), so subjects differ in how much and
#'               by which process they are contaminated. This is the only
#'               design row in which contamination varies between subjects,
#'               i.e. the only one that stresses the individual-differences
#'               estimands and the moderation analyses.
#' Methods (as user-facing pipelines): nothing (simple); hard cutoffs only,
#'          conventional (RT < 180 ms or > 3 s, the Wagenmakers et al., 2008,
#'          criteria) and aggressive (250 ms / 2.5 s); bmm only, as mixture
#'          down-weighting (ex-Gaussian / lognormal / inverse-Gaussian core) or
#'          as flag-and-exclude (threshold .5, or a Bernoulli draw with the
#'          probability as exclusion rate); combined pipelines, conventional
#'          cutoffs followed by mixture down-weighting or by flag-and-exclude;
#'          the ex-Gaussian mixture statistics with the accuracy channel
#'          corrected by adjust_ezdm_accuracy() (mixture_exgaussian_adjust;
#'          reuses the cell's fitted contaminant proportion, no extra EM);
#'          robust (median/IQR) and an oracle with perfect exclusion as
#'          references; and the three mixture methods reconstructed at the
#'          lower EM iteration caps 100 and 250.
#' Measures beyond bias/RMSE/correlations (all new in v7):
#'          excluded-set agreement between the threshold and the probabilistic
#'          exclusion rule (share of coinciding keep/drop decisions, plus the
#'          Jaccard overlap of the excluded sets); residual-location profile
#'          per exclusion pipeline (rate of surviving true contaminants,
#'          total, by generating type, and in the sub-ndt support-violation
#'          zone); eCDF calibration per exclusion pipeline (mean |delta eCDF|
#'          of the surviving trials against the same cell's pre-injection
#'          clean trials, plus the location-resolved deviation for an
#'          SBC-style band figure); mixed-row moderation (partial correlation
#'          of subject drift error with the zero-drift contaminant count,
#'          controlling the total count); and leading-edge depth moderation
#'          (correlation of subject-level flag sensitivity and drift error
#'          with the anticipation depth x_i).
#' Diagnostics: the EM is fitted once per core distribution per cell through
#'          flag_contaminant_rts(), whose "diagnostics" attribute carries the
#'          convergence flag and the iteration count that ezdm_summary_stats()
#'          discards. Because the EM is deterministic and maxit is only the loop
#'          bound, one fit at maxit = 500 yields the convergence rate at every
#'          lower cap, and - since bmm substitutes robust moments on
#'          non-convergence - the estimates at those caps as well. The same pass
#'          scores the flags against the generating contaminant labels
#'          (sensitivity, specificity, precision, F1) and supplies the
#'          probabilities the two exclusion rules use.
#' Runtime: ~4 h on 12 cores (20 rows vs v6's 22; the EM cost per cell is
#'          unchanged, and the ex-Gaussian EM dominates: its cost is bound by
#'          the number of EM iterations, not by the trial count, and
#'          per-process speed roughly halves under full parallel load). Use
#'          SIM_CORES to cap the workers on a loaded machine.

###############################################################################!
# 0) R Setup -------------------------------------------------------------------
###############################################################################!

pacman::p_load(here, bmm, tidyverse, SimDesign)
source(here("functions", "clean_plot.R"))

###############################################################################!
# 1) Simulation Design ---------------------------------------------------------
###############################################################################!

# The clean baseline (p_contam = 0) does not depend on the process, so it is
# run once per trial count under process = "none". The pure processes run at
# the two-point rate grid .02/.10 (mechanism diagnostics; the two points
# preserve the detection-collapse demonstration at reduced cost). The mixed
# process — the primary row — sweeps all three mean rates .02/.05/.10, with
# p_contam read as the *mean* of the subject-level Beta rate distribution.
# 2 + 3 x 2 x 2 + 3 x 2 = 20 rows.
Design <- createDesign(
  process  = c("none", "leading_edge", "delay", "guess", "mixed"),
  p_contam = c(0, 0.02, 0.05, 0.10),
  n_trials = c(50, 200),
  subset   = ((process == "none") == (p_contam == 0)) &
    (process == "mixed" | p_contam != 0.05)
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
  leading_edge_range = c(0.1, 0.5), # anticipation depth x_i ~ U(0.1, 0.5);
                                    # RT ~ U((1 - x_i) ndt_i, (1 + x_i) ndt_i)
  delay_range        = c(0, 2),     # delayed start-up: RT + U(0, 2 s)

  # Mixed process: each subject draws a total contamination rate
  # p_i ~ Beta(mean = p_contam, concentration = mixed_conc) and a process mix
  # w_i ~ Dirichlet(1, 1, 1), both shared across the subject's two conditions.
  # At mean .05 and concentration 10 the subject rates spread over roughly
  # 0-25%: many near-clean subjects plus a minority of heavy contaminators,
  # which is the leverage the mixed-row moderation analysis needs.
  mixed_conc = 10,

  # Hard pre-filtering cutoffs: the conventional criteria of the original
  # Wagenmakers et al. (2008) analysis (RT < 180 ms or > 3 s), plus an
  # aggressive variant (250 ms / 2.5 s) carried as a cutoffs-only method to
  # show the truncation trade-off.
  hard_lower = 0.18,
  hard_upper = 3.0,
  hard_aggr_lower = 0.25,
  hard_aggr_upper = 2.5,

  # Core distributions of the mixture method, and the EM iteration cap. bmm's
  # default (maxit = 100) is hit in a few percent of cells at tol = 1e-6, in
  # which case bmm falls back to robust moments with a warning; with 500 the
  # EM converges in > 99.8% of cells (the remaining fallbacks are counted).
  distributions = c("exgaussian", "lognormal", "invgaussian"),
  em_maxit      = 500,

  # Iteration caps at which convergence and estimates are scored. Only one EM
  # per core distribution is run, at em_maxit: .fit_rt_mixture() is
  # deterministic (method-of-moments initialisation, no RNG) and maxit is only
  # the loop bound `for (iter in seq_len(maxit))`, so a run capped at k <=
  # em_maxit would have converged iff `converged & iterations <= k`, and would
  # otherwise have returned bmm's robust-moment fallback. Both are recovered
  # from the single fit; see sections 4.2 and 4.5.
  em_maxit_grid = c(100, 250),

  # Contaminant-probability threshold for the two trial-exclusion rules.
  flag_threshold = 0.5,

  # Evaluation grid for the eCDF calibration: the surviving trials' eCDF is
  # read off at these quantiles of the same cell's pre-injection clean trials,
  # giving a location-resolved deviation (for the SBC-style band figure) and
  # its mean absolute value as a scalar summary.
  ecdf_probs = seq(0.05, 0.95, by = 0.05),

  # How the exclusion rules treat a cell whose EM did not converge, and whose
  # contaminant probabilities are therefore NA: "keep" retains every trial
  # (probability coalesced to 0), so the rule degrades to "no exclusion" for
  # that cell. This isolates the comparison of the two rules from the separate
  # question of EM convergence, whose rate is reported as clf_na. "drop"
  # removes the whole cell, which is what filter(contam_prob <= 0.5) does
  # silently, and is the defect this simulation documents.
  na_prob_rule = "keep"
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

# Correlations for the moderation measures. cor(use = "complete.obs") errors
# (rather than warns) when no complete pair exists — e.g. the subject-level
# flag sensitivity is NaN for every subject in the clean rows — so these
# return NA in that case and on zero variance, and NaN/Inf propagate through
# the partial correlation where the controls are collinear (pure rows, where
# one type's count equals the total count; allow_nan covers both).
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}

safe_pcor <- function(x, y, z) {
  rxy <- safe_cor(x, y)
  rxz <- safe_cor(x, z)
  ryz <- safe_cor(y, z)
  (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
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

  # Contamination profile per subject. In the pure-process rows every subject
  # has the same rate and a single process; in the mixed row each subject draws
  # a total rate p_i ~ Beta(mean = p_contam, concentration = mixed_conc) and a
  # process mix w_i ~ Dirichlet(1, 1, 1) (normalised Exp(1) draws), both shared
  # across the subject's two conditions, so contamination differs between
  # subjects in amount and in kind. The anticipation depth x_i ~ U(0.1, 0.5)
  # is drawn for every subject in every row (shared across the subject's two
  # conditions; unused where no leading-edge contaminants occur), so the
  # leading-edge mechanism reads the same x_i in the pure and the mixed rows.
  types <- c("leading_edge", "delay", "guess")
  x_i <- runif(fo$n_subjects, fo$leading_edge_range[1],
               fo$leading_edge_range[2])
  if (condition$process == "mixed") {
    w <- matrix(rexp(fo$n_subjects * 3), ncol = 3, dimnames = list(NULL, types))
    w <- w / rowSums(w)
    profile <- tibble(
      id  = 1:fo$n_subjects,
      p_i = rbeta(fo$n_subjects,
                  condition$p_contam * fo$mixed_conc,
                  (1 - condition$p_contam) * fo$mixed_conc),
      x_i = x_i
    ) |>
      bind_cols(as_tibble(w))
  } else {
    profile <- tibble(id = 1:fo$n_subjects, p_i = condition$p_contam,
                      x_i = x_i,
                      leading_edge = as.numeric(condition$process == "leading_edge"),
                      delay        = as.numeric(condition$process == "delay"),
                      guess        = as.numeric(condition$process == "guess"))
  }

  # Inject contaminants: round(n * p_i) trials per subject x condition cell
  data <- data_clean |>
    left_join(subj_params |> select(id, condition, bound, ndt),
              by = c("id", "condition")) |>
    left_join(profile, by = "id") |>
    group_by(id, condition) |>
    mutate(is_contaminant = sample(n()) <= round(n_tr * first(p_i))) |>
    ungroup()

  # Assign each contaminant a generating process (its subject's single process
  # in the pure rows; drawn from the subject's mix in the mixed row), then
  # apply the mechanisms per type. In pure rows only one type is non-empty.
  data$contam_type <- NA_character_
  idx <- which(data$is_contaminant)
  if (length(idx) > 0) {
    if (condition$process == "mixed") {
      data$contam_type[idx] <- vapply(idx, function(i) {
        sample(types, 1, prob = c(data$leading_edge[i], data$delay[i],
                                  data$guess[i]))
      }, character(1))
    } else {
      data$contam_type[idx] <- condition$process
    }

    idx_le <- idx[data$contam_type[idx] == "leading_edge"]
    if (length(idx_le) > 0) {
      data$rt[idx_le]       <- runif(length(idx_le),
                                     (1 - data$x_i[idx_le]) * data$ndt[idx_le],
                                     (1 + data$x_i[idx_le]) * data$ndt[idx_le])
      data$response[idx_le] <- rbinom(length(idx_le), 1, 0.5)
    }
    idx_dl <- idx[data$contam_type[idx] == "delay"]
    if (length(idx_dl) > 0) {
      data$rt[idx_dl] <- data$rt[idx_dl] +
        runif(length(idx_dl), fo$delay_range[1], fo$delay_range[2])
    }
    idx_gu <- idx[data$contam_type[idx] == "guess"]
    if (length(idx_gu) > 0) {
      redraw <- rddm(length(idx_gu), drift = 0, bound = data$bound[idx_gu],
                     ndt = data$ndt[idx_gu], zr = 0.5)
      data$rt[idx_gu]       <- redraw$rt
      data$response[idx_gu] <- redraw$response
    }
  }

  # Truth carried per subject x condition: generating parameters, the
  # population moments of the core process, the subject's true contamination
  # rate (constant within pure rows), the process mix shares (indicator
  # columns in pure rows), and the anticipation depth x_i
  truth <- subj_params |>
    rename(drift_true = drift, bound_true = bound, ndt_true = ndt) |>
    mutate(ez_forward(drift_true, bound_true, ndt_true)) |>
    left_join(profile |> select(id, p_contam_true = p_i, x_i,
                                w_leading_edge = leading_edge,
                                w_delay = delay, w_guess = guess),
              by = "id")

  # Contaminants *replace* trials, so the pre-injection clean versions are
  # otherwise lost; they are returned as the eCDF calibration target.
  list(
    data  = data |> select(id, condition, trial, rt, response,
                           is_contaminant, contam_type),
    clean = data_clean,
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
  clean <- dat$clean
  truth <- dat$truth
  pop   <- dat$pop

  # Hard pre-filtered versions: the conventional cutoffs (applied before any
  # method, and the basis of the combined pipelines) and the aggressive
  # cutoffs (carried as a cutoffs-only method)
  d_pf      <- d |> filter(rt >= fo$hard_lower, rt <= fo$hard_upper)
  d_pf_aggr <- d |> filter(rt >= fo$hard_aggr_lower, rt <= fo$hard_aggr_upper)

  ## 4.1) EM pass: one fit per core distribution per subject x condition cell --
  # ezdm_summary_stats() discards the EM's convergence flag and iteration
  # count, so they are harvested from flag_contaminant_rts(), which attaches
  # them as a "diagnostics" attribute. For version = "3par" the two functions
  # fit the *same* EM: both run it on the pooled rt vector of the cell with
  # bounds resolved from that same vector (bmm helpers-data.R, l. 466/493 and
  # l. 803), and both leave init_contaminant, max_contaminant and tol at their
  # defaults. The iteration counts therefore describe the mixture fits used
  # below, and the ex-Gaussian probabilities are the ones a user filtering on
  # contam_prob would obtain.
  #
  # The three cores are timed back to back inside the same cell, so they share
  # the same worker and the same instantaneous parallel load; only the ratio
  # between them is interpretable, never the absolute seconds. Mean iterations
  # is carried alongside as a load-independent cost proxy.
  em_cell <- function(cell) {
    fits <- map(fo$distributions, function(dist) {
      t0 <- proc.time()
      p  <- suppressWarnings(flag_contaminant_rts(
        cell$rt, distribution = dist, contaminant_bound = c("min", "max"),
        maxit = fo$em_maxit
      ))
      dt <- proc.time() - t0
      dg <- attr(p, "diagnostics")
      list(prob = as.numeric(p),
           diag = tibble(id = cell$id[1], condition = cell$condition[1],
                         distribution = dist,
                         converged  = dg$converged,
                         iterations = dg$iterations,
                         em_secs = dt[["elapsed"]], em_cpu = dt[["user.self"]]))
    })
    names(fits) <- fo$distributions
    list(diag  = map_dfr(fits, "diag"),
         trial = cell |> mutate(contam_prob = fits[["exgaussian"]]$prob))
  }

  # group_split() keeps cell$rt and the returned probability vector aligned
  # element-wise, which summarize() would not guarantee.
  em      <- d |> group_by(id, condition) |> group_split() |> map(em_cell)
  em_diag <- map_dfr(em, "diag")
  d_prob  <- map_dfr(em, "trial")

  ## 4.2) Trial sets for the two exclusion rules -------------------------------
  # Both reuse the ex-Gaussian probabilities from 4.1; no extra EM is run.
  p_ex <- if (fo$na_prob_rule == "keep") coalesce(d_prob$contam_prob, 0) else
    replace_na(d_prob$contam_prob, 1)

  drop_thr <- p_ex > fo$flag_threshold
  d_thr    <- d_prob[!drop_thr, ]

  ## 4.2b) Combined pipeline: hard cutoffs, then flag and exclude --------------
  # Mechanistic rationale: the EM's sensitivity to fast guesses collapses at
  # high rates because the fitted core stretches to absorb them; removing the
  # sub-cutoff spike before the EM sees it should let the core fit cleanly and
  # restore detection of what remains. This requires one extra ex-Gaussian EM
  # per cell, on the pre-filtered data.
  flag_cell_pf <- function(cell) {
    p <- suppressWarnings(flag_contaminant_rts(
      cell$rt, distribution = "exgaussian",
      contaminant_bound = c("min", "max"), maxit = fo$em_maxit
    ))
    cell |> mutate(contam_prob = as.numeric(p))
  }
  d_pf_prob <- d_pf |> group_by(id, condition) |> group_split() |>
    map_dfr(flag_cell_pf)
  p_ex_pf <- if (fo$na_prob_rule == "keep")
    coalesce(d_pf_prob$contam_prob, 0) else replace_na(d_pf_prob$contam_prob, 1)
  d_thr_pf <- d_pf_prob[p_ex_pf <= fo$flag_threshold, ]

  # SimDesign does not re-seed between replications: mainsim() reads
  # .GlobalEnv$.Random.seed and lets the worker's stream run on, so RNG drawn
  # here would shift the data Generate() produces in every later replication
  # handled by this worker. The stream is saved here and restored after the
  # last stochastic analysis step — the accuracy-adjustment draws in 4.4; the
  # summary-statistic methods in between are deterministic — so that
  # Generate() sees exactly the state it saw before these sections existed.
  seed_before <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  drop_prb <- rbinom(nrow(d_prob), 1, p_ex) == 1
  d_prb    <- d_prob[!drop_prb, ]

  ## 4.2c) Excluded-set agreement between the two exclusion rules --------------
  # Do the threshold rule and the probabilistic rule remove the same trials?
  # agr_share is the per-condition share of trials on which the keep/drop
  # decisions coincide; agr_jaccard is the overlap among the excluded trials
  # (|both| / |either|; NaN when neither rule excludes anything, as in the
  # clean rows). This is what licenses — or refuses — "interchangeable" for
  # any downstream model that consumes the surviving trials.
  agr <- tibble(condition = d_prob$condition,
                thr = drop_thr, prb = drop_prb) |>
    summarize(
      agr_share   = mean(thr == prb),
      agr_jaccard = sum(thr & prb) / sum(thr | prb),
      .by = condition
    ) |>
    mutate(method = "excl_agreement")

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
    sumstats(d_pf, "mixture", "mixture_exgaussian_pf", "exgaussian"),
    sumstats(d_thr, "simple", "flag_threshold"),
    sumstats(d_prb, "simple", "flag_probabilistic"),
    sumstats(d_pf_aggr, "simple", "simple_pf_aggr"),
    sumstats(d_thr_pf, "simple", "flag_threshold_pf")
  )

  ## 4.4) Accuracy-adjusted mixture statistics ---------------------------------
  # mixture_exgaussian_adjust = the ex-Gaussian mixture statistics with
  # adjust_ezdm_accuracy() applied: the response counts are corrected for the
  # estimated share of chance-accuracy contaminants, reusing the cell's fitted
  # contaminant proportion (no extra EM; guess_rate stays at the accuracy-
  # coding chance level 0.5). The function is stochastic — it draws the
  # contaminant count Binomial(n_trials, pi_c) and the chance-correct count
  # Binomial(n_contam, 0.5) — which is why it sits inside the seed guard, and
  # scalar-only, so it is applied per cell. A cell whose EM did not converge
  # (contaminant_prop NA) is left unchanged by the function's own NA rule.
  # The RT moments are the mixture's; only n_upper and n_trials change, so
  # the closed-form inversion below sees the corrected proportion correct.
  mix_ref <- ss |> filter(method == "mixture_exgaussian")
  adj <- pmap(list(mix_ref$n_upper, mix_ref$n_trials, mix_ref$contaminant_prop),
              \(nu, nt, pc) adjust_ezdm_accuracy(nu, nt, pc)) |>
    bind_rows()
  ss <- bind_rows(
    ss,
    mix_ref |> mutate(n_upper  = adj$n_upper_adj,
                      n_trials = adj$n_trials_adj,
                      method   = "mixture_exgaussian_adjust")
  )

  if (!is.null(seed_before)) assign(".Random.seed", seed_before, envir = .GlobalEnv)

  ## 4.5) Reconstructed mixture statistics at lower iteration caps -------------
  # What ezdm_summary_stats(maxit = k) would have returned, for k < em_maxit.
  # bmm substitutes the robust moments whenever the EM does not converge
  # (helpers-data.R, l. 481-484), and the EM is deterministic, so a cell's
  # result at cap k is its mixture result if it converged within k iterations
  # and its robust result otherwise. Both are already in `ss`; nothing is
  # refitted. Verified against ezdm_summary_stats(maxit = 100) on both branches.
  rob_ref <- ss |>
    filter(method == "robust") |>
    select(id, condition, mean_rt_rob = mean_rt, var_rt_rob = var_rt)

  ss_capped <- map_dfr(fo$em_maxit_grid, function(k) {
    map_dfr(fo$distributions, function(dist) {
      ss |>
        filter(method == paste0("mixture_", dist)) |>
        left_join(em_diag |> filter(distribution == dist) |>
                    select(id, condition, converged, iterations),
                  by = c("id", "condition")) |>
        left_join(rob_ref, by = c("id", "condition")) |>
        mutate(
          kept    = converged & iterations <= k,
          mean_rt = if_else(kept, mean_rt, mean_rt_rob),
          var_rt  = if_else(kept, var_rt, var_rt_rob),
          contaminant_prop = if_else(kept, contaminant_prop, NA_real_),
          method  = paste0("mixture_", dist, "_it", k)
        ) |>
        select(-converged, -iterations, -kept, -mean_rt_rob, -var_rt_rob)
    })
  })

  ss <- bind_rows(ss, ss_capped)

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

  # Per-cell contaminant counts for the mixed-row moderation: how many of the
  # cell's trials were contaminated in total, and by the zero-drift process
  contam_counts <- d |>
    summarize(
      n_contam_i = sum(is_contaminant),
      n_guess_i  = sum(is_contaminant & contam_type == "guess"),
      .by = c(id, condition)
    )

  per_cell <- est |>
    left_join(oracle_ref, by = c("id", "condition")) |>
    left_join(contam_counts, by = c("id", "condition")) |>
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
      # Differential-contamination damage: how strongly a subject's drift error
      # tracks their true contamination rate. Defined only in the mixed row,
      # where rates vary between subjects (NaN elsewhere, where sd(p) = 0).
      cor_pcon  = suppressWarnings(cor(e_drift, p_contam_true,
                                       use = "complete.obs")),
      # Mixed-row moderation (pre-committed, spec 2026-08-21 §3.4): partial
      # correlation of the subject's drift error with their zero-drift
      # contaminant count, controlling the total contaminant count. If it is
      # substantial (> 2 x MCSE), the residual differential bias is carried
      # by the undetectable share specifically, which scopes the
      # differential-bias claim to undetectable contamination. NaN outside
      # the mixed row (counts collinear or constant).
      mod_pcor_guess = safe_pcor(e_drift, n_guess_i, n_contam_i),
      # Leading-edge depth moderation: does the subject's drift error track
      # the anticipation depth x_i? Interpretable in the leading_edge rows;
      # computed everywhere so every replication returns the same key set.
      dep_cor_edrift = safe_cor(e_drift, x_i),
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

  ## 4.6) EM convergence and cost by core distribution -------------------------
  # `converged & iterations <= k` reproduces exactly what a run capped at k
  # would have reported. em_break isolates a numerical early break (NA
  # responsibilities or a non-finite log-likelihood, which stop the EM without
  # setting converged) from em_capped, i.e. genuinely running out of iterations.
  em_by <- em_diag |>
    summarize(
      em_conv100   = mean(converged & iterations <= 100),
      em_conv250   = mean(converged & iterations <= 250),
      em_conv500   = mean(converged & iterations <= fo$em_maxit),
      em_iter_mean = mean(iterations),
      em_iter_med  = median(iterations),
      em_break     = mean(!converged & iterations <  fo$em_maxit),
      em_capped    = mean(!converged & iterations >= fo$em_maxit),
      em_secs      = mean(em_secs),
      em_cpu       = mean(em_cpu),
      .by = c(distribution, condition)
    ) |>
    mutate(method = paste0("em_", distribution), .keep = "unused")

  ## 4.7) Flagging performance against ground truth ----------------------------
  # Counts are pooled over subjects within a condition before the rates are
  # formed. A per-subject rate is not usable here: at 50 trials and 2%
  # contamination a cell holds a single contaminant, so its sensitivity is 0
  # or 1, and it is undefined on clean data. The per-process breakdown comes
  # from the design factor `process`; the p_contam = 0 rows give the
  # false-positive rate on wholly clean data (1 - clf_spec).
  clf <- d_prob |>
    mutate(flagged = contam_prob > fo$flag_threshold) |>
    summarize(
      clf_tp = sum(is_contaminant  &  flagged, na.rm = TRUE),
      clf_fp = sum(!is_contaminant &  flagged, na.rm = TRUE),
      clf_tn = sum(!is_contaminant & !flagged, na.rm = TRUE),
      clf_fn = sum(is_contaminant  & !flagged, na.rm = TRUE),
      clf_na = mean(is.na(contam_prob)),
      # Type-resolved sensitivity: in the pure rows only the row's own type is
      # defined (the others are 0/0 = NaN); in the mixed row all three are,
      # which shows which components of the cocktail the flags actually reach.
      clf_sens_le = sum(flagged & contam_type == "leading_edge", na.rm = TRUE) /
        sum(contam_type == "leading_edge", na.rm = TRUE),
      clf_sens_dl = sum(flagged & contam_type == "delay", na.rm = TRUE) /
        sum(contam_type == "delay", na.rm = TRUE),
      clf_sens_gu = sum(flagged & contam_type == "guess", na.rm = TRUE) /
        sum(contam_type == "guess", na.rm = TRUE),
      .by = condition
    ) |>
    mutate(
      clf_sens = clf_tp / (clf_tp + clf_fn),
      clf_spec = clf_tn / (clf_tn + clf_fp),
      clf_prec = clf_tp / (clf_tp + clf_fp),
      clf_f1   = 2 * clf_prec * clf_sens / (clf_prec + clf_sens),
      method   = "flag_exgaussian"
    )

  # Depth moderation of the flags (leading_edge rows): per subject, the share
  # of their true contaminants flagged, correlated with the anticipation
  # depth x_i. Deep anticipations (large x_i) reach further into the core, so
  # a negative correlation means the flags miss exactly the subjects whose
  # estimates are damaged most. sens_i is NaN for subjects without
  # contaminants (all of them in the clean rows), which safe_cor() drops.
  subj_sens <- d_prob |>
    mutate(flagged = contam_prob > fo$flag_threshold) |>
    summarize(
      sens_i = sum(flagged & is_contaminant, na.rm = TRUE) /
        sum(is_contaminant),
      .by = c(id, condition)
    ) |>
    left_join(truth |> select(id, condition, x_i), by = c("id", "condition"))

  clf <- clf |>
    left_join(
      subj_sens |>
        summarize(dep_cor_sens = safe_cor(sens_i, x_i), .by = condition),
      by = "condition"
    )

  ## 4.8) Residual-location profile per exclusion pipeline ---------------------
  # What the trial-level reader's risk actually is: among the trials a
  # pipeline keeps, the rate of true contaminants — in total, by generating
  # type, and in the sub-ndt support-violation zone (rt < ndt_i), where a
  # surviving trial has zero Wiener likelihood and forces ndt below the
  # fastest surviving RT, rather than merely leaning on a tail. The cutoff
  # pipelines are the surviving sets of simple_pf / simple_pf_aggr, so the
  # measures attach to those method names.
  pipes <- list(simple_pf          = d_pf,
                simple_pf_aggr     = d_pf_aggr,
                flag_threshold     = d_thr,
                flag_probabilistic = d_prb,
                flag_threshold_pf  = d_thr_pf)

  ndt_ref <- truth |> select(id, condition, ndt_true)

  resloc <- imap(pipes, \(dd, label) {
    dd |>
      left_join(ndt_ref, by = c("id", "condition")) |>
      summarize(
        res_contam = mean(is_contaminant),
        res_subndt = mean(is_contaminant & rt < ndt_true),
        res_le = mean(is_contaminant & contam_type == "leading_edge"),
        res_dl = mean(is_contaminant & contam_type == "delay"),
        res_gu = mean(is_contaminant & contam_type == "guess"),
        .by = condition
      ) |>
      mutate(method = label)
  }) |>
    list_rbind()

  ## 4.9) eCDF calibration per exclusion pipeline ------------------------------
  # Distributional calibration of the surviving trials against the same
  # cell's pre-injection clean trials: the surviving eCDF is read off at the
  # clean cell's quantiles (fo$ecdf_probs), giving a location-resolved signed
  # deviation (positive = surplus mass below that location) whose mean over
  # cells feeds the SBC-style band figure, plus the mean |deviation| per cell
  # as the scalar summary (ecdf_mad).
  clean_cells <- clean |>
    summarize(rt_clean = list(rt), .by = c(id, condition))

  ecdf_dev_cell <- function(rt_surv, rt_clean) {
    q <- quantile(rt_clean, probs = fo$ecdf_probs, names = FALSE)
    if (length(rt_surv) == 0) return(rep(NA_real_, length(q)))
    ecdf(rt_surv)(q) - fo$ecdf_probs
  }

  ecdf_cal <- imap(pipes, \(dd, label) {
    cells <- dd |>
      summarize(rt_surv = list(rt), .by = c(id, condition)) |>
      left_join(clean_cells, by = c("id", "condition")) |>
      mutate(dev = map2(rt_surv, rt_clean, ecdf_dev_cell))
    bind_rows(
      cells |>
        mutate(value = map_dbl(dev, \(x) mean(abs(x)))) |>
        summarize(value = mean(value, na.rm = TRUE), .by = condition) |>
        mutate(measure = "ecdf_mad"),
      cells |>
        mutate(dev = map(dev, \(x)
          tibble(measure = sprintf("ecdf_dev_p%02d",
                                   round(100 * fo$ecdf_probs)),
                 value = x))) |>
        select(condition, dev) |>
        unnest(dev) |>
        summarize(value = mean(value, na.rm = TRUE),
                  .by = c(condition, measure))
    ) |>
      mutate(method = label)
  }) |>
    list_rbind()

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
      pivot_longer(c(bias_effect, mse_effect), names_to = "measure"),
    em_by |>
      pivot_longer(-c(method, condition), names_to = "measure"),
    clf |>
      pivot_longer(-c(method, condition), names_to = "measure"),
    agr |>
      pivot_longer(-c(method, condition), names_to = "measure"),
    resloc |>
      pivot_longer(-c(method, condition), names_to = "measure"),
    ecdf_cal
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

  # Monte Carlo SE for every bias measure; for the flagging and EM
  # convergence rates, which the appendix reports as point estimates; and for
  # the v7 measures — exclusion-rule agreement (agr_), residual-location
  # profile (res_), eCDF calibration (ecdf_, including every location on the
  # grid: the SBC-style band figure is mean +/- 2 MCSE), and the moderation
  # correlations (mod_, dep_), which the pre-committed decision rules read
  # against 2 x MCSE.
  is_mcse <- str_detect(colnames(results),
                        "^(bias_|clf_(sens|spec|prec|f1)|em_conv|agr_|res_|ecdf_|mod_|dep_)")
  R    <- results[, is_mcse, drop = FALSE]
  mcse <- apply(R, 2, sd, na.rm = TRUE) / sqrt(colSums(!is.na(R)))
  names(mcse) <- if_else(str_starts(names(mcse), "bias_"),
                         str_replace(names(mcse), "^bias_", "mcse_"),
                         paste0("mcse_", names(mcse)))

  c(m, mcse)
}

###############################################################################!
# 6) Run Simulation ------------------------------------------------------------
###############################################################################!

n_reps  <- 200

# Leave headroom for anything else running on the machine (model fitting in a
# second R session, in particular): set SIM_CORES to cap the workers, e.g.
# SIM_CORES=8 Rscript scripts/02_simulation_contaminant_methods.R
n_cores <- min(parallel::detectCores() - 1, 15)
if (nzchar(Sys.getenv("SIM_CORES"))) {
  n_cores <- min(n_cores, as.integer(Sys.getenv("SIM_CORES")))
}

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
                   "mixture_lognormal", "mixture_invgaussian",
                   "mixture_exgaussian_adjust", "oracle",
                   "simple_pf", "robust_pf", "mixture_exgaussian_pf",
                   # Trial-exclusion rules and the diagnostic pseudo-methods
                   "flag_threshold", "flag_probabilistic",
                   "simple_pf_aggr", "flag_threshold_pf",
                   "mixture_exgaussian_it100", "mixture_lognormal_it100",
                   "mixture_invgaussian_it100",
                   "mixture_exgaussian_it250", "mixture_lognormal_it250",
                   "mixture_invgaussian_it250",
                   "em_exgaussian", "em_lognormal", "em_invgaussian",
                   "flag_exgaussian", "excl_agreement")
method_labels <- c("Simple", "Robust", "Mixture (ex-Gaussian)",
                   "Mixture (lognormal)", "Mixture (inv-Gaussian)",
                   "Mixture (ex-Gaussian) + accuracy adjustment",
                   "Oracle (perfect exclusion)",
                   "Simple, pre-filtered", "Robust, pre-filtered",
                   "Mixture (ex-Gaussian), pre-filtered",
                   "Threshold exclusion (p > .5)", "Probabilistic exclusion",
                   "Simple, aggressive pre-filter",
                   "Pre-filter + threshold exclusion",
                   "Mixture (ex-Gaussian), maxit 100",
                   "Mixture (lognormal), maxit 100",
                   "Mixture (inv-Gaussian), maxit 100",
                   "Mixture (ex-Gaussian), maxit 250",
                   "Mixture (lognormal), maxit 250",
                   "Mixture (inv-Gaussian), maxit 250",
                   "EM (ex-Gaussian)", "EM (lognormal)", "EM (inv-Gaussian)",
                   "Flagging (ex-Gaussian)", "Exclusion-rule agreement")
process_labels <- c(leading_edge = "Leading-edge anticipations (around ndt)",
                    delay = "Delayed start-ups (+ U(0, 2 s))",
                    guess = "Informationless responses (drift 0)",
                    mixed = "Mixed processes (subject-specific rates)")

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

results_long <- bind_rows(results_long |> filter(process != "none"), baseline)

# factor() silently converts an unregistered method to NA, which would drop it
# from the CSV without any error, so check the mapping is complete first.
stopifnot(all(unique(results_long$method) %in% method_levels))

results_long <- results_long |>
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
                  "Mixture (ex-Gaussian) + accuracy adjustment",
                  "Oracle (perfect exclusion)")
main_colors <- c(
  "Simple" = "grey50",
  "Robust" = eam_palette[2],
  "Mixture (ex-Gaussian)" = eam_palette[1],
  "Mixture (ex-Gaussian) + accuracy adjustment" = eam_palette[3],
  "Oracle (perfect exclusion)" = "black"
)

# The rate-gradient figures show the three pure processes at the two-point
# rate grid (plus the shared clean baseline); the mixed results get their own
# pipeline figure (8.3) with the mean-rate sweep as a facet dimension.
pure_processes <- process_labels[c("leading_edge", "delay", "guess")]

fig1_data <- results_long |>
  filter(n_trials == 200, method %in% main_methods, process %in% pure_processes,
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
  scale_x_continuous(labels = scales::percent, breaks = c(0, .02, .10)) +
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
  filter(n_trials == 200, method %in% core_methods, process %in% pure_processes,
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
  scale_x_continuous(labels = scales::percent, breaks = c(0, .02, .10)) +
  labs(x = "Contamination rate", y = "Bias (estimate − generating value)",
       color = "Core distribution", linetype = "Condition") +
  clean_plot() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 8), legend.spacing.y = unit(0, "pt"))

ggsave(here("figures", "02_mixture_core_distributions.pdf"),
       p_fig2, width = 6.5, height = 6, device = "pdf")

## 8.3) Figure 3: preprocessing pipelines under mixed contamination ------------
# The primary comparison: nothing / cutoffs only / bmm only / combined, on
# the realistic cocktail with subject-specific rates and mixes, across the
# mean-rate sweep (facet rows).
pipeline_methods <- c(
  "Simple", "Simple, pre-filtered", "Simple, aggressive pre-filter",
  "Robust", "Mixture (ex-Gaussian)",
  "Mixture (ex-Gaussian) + accuracy adjustment",
  "Mixture (ex-Gaussian), pre-filtered",
  "Threshold exclusion (p > .5)", "Pre-filter + threshold exclusion",
  "Oracle (perfect exclusion)"
)

rate_label <- \(x) paste0("Mean rate ", scales::percent(as.numeric(x)))

fig3_data <- results_long |>
  filter(process == process_labels["mixed"], p_contam > 0,
         method %in% pipeline_methods,
         measure %in% c("bias_drift", "bias_bound", "bias_ndt")) |>
  mutate(parameter = factor(str_remove(measure, "bias_"),
                            levels = names(param_labels),
                            labels = param_labels),
         method = factor(method, levels = rev(pipeline_methods)))

p_fig3 <- fig3_data |>
  ggplot(aes(x = value, y = method, color = condition,
             shape = factor(n_trials))) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  facet_grid(p_contam ~ parameter, scales = "free_x",
             labeller = labeller(p_contam = rate_label)) +
  scale_color_manual(values = c("Easy" = eam_palette[1],
                                "Hard" = eam_palette[2])) +
  labs(x = "Bias (estimate − generating value)", y = NULL,
       color = "Condition", shape = "Trials per cell") +
  clean_plot() +
  theme(legend.position = "bottom")

ggsave(here("figures", "02_mixed_pipelines.pdf"),
       p_fig3, width = 6.5, height = 8.5, device = "pdf")

## 8.4) Figure 4: eCDF calibration, SBC-style (mixed process, n = 200) ---------
# Laid out like the standard simulation-based-calibration ECDF-difference plot
# (Talts et al., 2018; Sailynoja, Burkner & Vehtari, 2022): the horizontal axis
# is the normalised rank statistic z, the vertical axis the deviation of the
# empirical CDF from uniformity, F(z) - z, and a shaded envelope marks the
# region in which a calibrated curve is expected to stay.
#
# The plotted quantity is exactly F(z) - z: the surviving trials' eCDF read at
# the z-quantile of the same cell's pre-injection clean trials, minus z. A
# pipeline that left the surviving RT distribution indistinguishable from clean
# would sit on the dashed zero line.
#
# Two envelopes, which must not be confused:
#   * the pale background lens is the *single-cell reference*, 1.96 sqrt(z(1-z)/n)
#     for n trials in one subject x condition cell. It is the sampling noise a
#     reader's own single cell would show under perfect calibration, and it is
#     the envelope the traditional SBC plot draws.
#   * the narrow coloured ribbon is +/- 2 MCSE of the plotted mean, i.e. the
#     Monte Carlo precision of this simulation. It is far tighter than the lens
#     because the mean is taken over cells and replications, so a curve can sit
#     well inside the lens and still be reliably non-zero.
exclusion_methods <- c("Simple, pre-filtered", "Simple, aggressive pre-filter",
                       "Threshold exclusion (p > .5)",
                       "Probabilistic exclusion",
                       "Pre-filter + threshold exclusion")

fig4_data <- results_long |>
  filter(process == process_labels["mixed"], n_trials == 200,
         method %in% exclusion_methods,
         str_detect(measure, "^(mcse_)?ecdf_dev_p")) |>
  mutate(stat = if_else(str_starts(measure, "mcse_"), "mcse", "dev"),
         location = as.numeric(str_extract(measure, "\\d+$")) / 100) |>
  select(-measure) |>
  pivot_wider(names_from = stat, values_from = value) |>
  mutate(method = factor(method, levels = exclusion_methods))

# Single-cell reference lens, drawn on a fine grid so the ends close smoothly
sbc_lens <- tibble(location = seq(0, 1, length.out = 401)) |>
  mutate(half = 1.96 * sqrt(location * (1 - location) / 200))

fig4_colors <- setNames(eam_palette[1:length(exclusion_methods)],
                        exclusion_methods)

p_fig4 <- ggplot() +
  geom_ribbon(data = sbc_lens,
              aes(x = location, ymin = -half, ymax = half),
              fill = "#DCE8F7", color = "#9DB8D9", linewidth = 0.25) +
  geom_hline(yintercept = 0, color = "grey35", linetype = "dashed",
             linewidth = 0.4) +
  geom_ribbon(data = fig4_data,
              aes(x = location, ymin = dev - 2 * mcse, ymax = dev + 2 * mcse,
                  fill = method),
              alpha = 0.3, color = NA) +
  geom_line(data = fig4_data,
            aes(x = location, y = dev, color = method), linewidth = 0.7) +
  facet_grid(p_contam ~ condition, labeller = labeller(p_contam = rate_label)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(0.02)) +
  scale_color_manual(values = fig4_colors) +
  scale_fill_manual(values = fig4_colors) +
  labs(x = "Normalized rank statistic (z)",
       y = "ECDF difference   F(z) − z",
       color = "Exclusion pipeline", fill = "Exclusion pipeline") +
  # Three rows (two columns): the pipeline names are too long to fit three
  # columns at 6.5 in, which truncates the last label.
  guides(color = guide_legend(nrow = 3, byrow = TRUE, title.position = "top"),
         fill  = guide_legend(nrow = 3, byrow = TRUE, title.position = "top")) +
  clean_plot() +
  theme(panel.border = element_rect(color = "grey25", fill = NA,
                                    linewidth = 0.4),
        panel.spacing.x = unit(12, "pt"),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(size = 8),
        legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 8),
        legend.key.width = unit(18, "pt"))

ggsave(here("figures", "02_ecdf_calibration.pdf"),
       p_fig4, width = 6.5, height = 7.5, device = "pdf")

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

cat("\n=== 10% contamination, n = 200: bias in EZ parameters by process ===\n\n")
results_long |>
  filter(p_contam == 0.10, n_trials == 200, condition != "All",
         measure %in% c("bias_drift", "bias_bound", "bias_ndt", "pi_c")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(process, method, condition, bias_drift, bias_bound, bias_ndt, pi_c) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== 10% contamination, n = 200: practical estimands ===\n\n")
results_long |>
  filter(p_contam == 0.10, n_trials == 200,
         measure %in% c("bias_effect", "rmse_effect", "cor_drift",
                        "cor_bound", "cor_ndt")) |>
  pivot_wider(names_from = c(measure, condition), values_from = value) |>
  select(process, method, bias_effect_All, rmse_effect_All,
         cor_drift_Easy, cor_drift_Hard, cor_bound_Easy, cor_ndt_Easy) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== Pre-filtering: bias vs the unfiltered method (10%, n = 200) ===\n\n")
results_long |>
  filter(p_contam == 0.10, n_trials == 200, condition != "All",
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

cat("\n=== EM convergence and cost by core distribution ===\n\n")
cat("em_secs / em_cpu are recorded under full parallel load and are only\n",
    "comparable between cores within a row; em_iter_mean is load-independent.\n\n")
results_long |>
  filter(str_starts(as.character(method), "EM "), condition != "All") |>
  summarize(value = mean(value), .by = c(measure, n_trials, method)) |>
  pivot_wider(names_from = measure, values_from = value) |>
  transmute(n_trials, method,
            conv_100 = round(em_conv100, 4), conv_250 = round(em_conv250, 4),
            conv_500 = round(em_conv500, 4),
            iter_mean = round(em_iter_mean, 1), iter_med = round(em_iter_med, 1),
            early_break = signif(em_break, 3), capped = signif(em_capped, 3),
            cpu_ms = round(1000 * em_cpu, 2),
            cpu_rel = round(em_cpu / min(em_cpu), 2)) |>
  arrange(n_trials, method) |>
  print(n = Inf)

cat("\n=== Effect of the iteration cap on drift and bound bias (10%, n = 200) ===\n\n")
results_long |>
  filter(p_contam == 0.10, n_trials == 200, condition != "All",
         str_starts(as.character(method), "Mixture"),
         !str_detect(as.character(method), "pre-filtered"),
         measure %in% c("bias_drift", "bias_bound", "fallback")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(process, method, condition, bias_drift, bias_bound, fallback) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== Flagging vs. ground truth (ex-Gaussian core, threshold .5) ===\n\n")
results_long |>
  filter(method == "Flagging (ex-Gaussian)", condition != "All",
         measure %in% c("clf_sens", "clf_spec", "clf_prec", "clf_f1", "clf_na")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(n_trials, process, p_contam, condition,
         clf_sens, clf_spec, clf_prec, clf_f1, clf_na) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(n_trials, process, p_contam, condition) |>
  print(n = Inf)

cat("\n=== Exclusion rules vs. down-weighting (10%, n = 200) ===\n\n")
results_long |>
  filter(p_contam == 0.10, n_trials == 200, condition != "All",
         measure %in% c("bias_drift", "bias_bound", "bias_ndt"),
         method %in% c("Simple", "Mixture (ex-Gaussian)",
                       "Threshold exclusion (p > .5)", "Probabilistic exclusion",
                       "Oracle (perfect exclusion)")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(process, method, condition, bias_drift, bias_bound, bias_ndt) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(n = Inf)

cat("\n=== Excluded-set agreement: threshold vs probabilistic rule ===\n\n")
results_long |>
  filter(method == "Exclusion-rule agreement", condition != "All",
         measure %in% c("agr_share", "agr_jaccard")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(n_trials, process, p_contam, condition, agr_share, agr_jaccard) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(n_trials, process, p_contam, condition) |>
  print(n = Inf)

cat("\n=== Residual-location profile: surviving true contaminants ===\n\n")
cat("res_contam = rate among surviving trials; res_subndt = of these, in the\n",
    "sub-ndt support-violation zone; res_le/dl/gu = by generating type.\n\n")
results_long |>
  filter(p_contam > 0, n_trials == 200, condition != "All",
         measure %in% c("res_contam", "res_subndt", "res_le",
                        "res_dl", "res_gu")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(process, p_contam, method, condition, res_contam, res_subndt,
         res_le, res_dl, res_gu) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(process, p_contam, method, condition) |>
  print(n = Inf)

cat("\n=== eCDF calibration: mean |deviation| vs the clean trials ===\n\n")
results_long |>
  filter(p_contam > 0, n_trials == 200, condition != "All",
         measure == "ecdf_mad") |>
  pivot_wider(names_from = condition, values_from = value) |>
  select(process, p_contam, method, Easy, Hard) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(process, p_contam, method) |>
  print(n = Inf)

cat("\n=== Accuracy adjustment vs plain mixture (adjust decision rule) ===\n\n")
results_long |>
  filter(p_contam > 0,
         method %in% c("Mixture (ex-Gaussian)",
                       "Mixture (ex-Gaussian) + accuracy adjustment"),
         measure %in% c("bias_drift", "mcse_drift", "bias_pc",
                        "bias_effect", "cor_drift")) |>
  pivot_wider(names_from = c(measure, condition), values_from = value) |>
  transmute(n_trials, process, p_contam, method,
            drift_easy = round(bias_drift_Easy, 3),
            mcse_easy  = round(mcse_drift_Easy, 4),
            drift_hard = round(bias_drift_Hard, 3),
            pc_easy    = round(bias_pc_Easy, 4),
            effect     = round(bias_effect_All, 3),
            cor_drift_easy = round(cor_drift_Easy, 3)) |>
  arrange(n_trials, process, p_contam, method) |>
  print(n = Inf)

cat("\n=== Mixed-row moderation: drift error vs zero-drift count ===\n\n")
cat("Partial correlation of subject drift error with the zero-drift\n",
    "contaminant count, controlling the total count (mean +/- MCSE);\n",
    "the pre-committed rule reads it against 2 x MCSE.\n\n")
results_long |>
  filter(process == process_labels["mixed"], p_contam > 0, condition != "All",
         measure %in% c("mod_pcor_guess", "mcse_mod_pcor_guess"),
         method %in% c("Simple", "Robust", "Mixture (ex-Gaussian)",
                       "Mixture (ex-Gaussian) + accuracy adjustment",
                       "Threshold exclusion (p > .5)",
                       "Pre-filter + threshold exclusion",
                       "Oracle (perfect exclusion)")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(n_trials, p_contam, method, condition,
         pcor = mod_pcor_guess, mcse = mcse_mod_pcor_guess) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(n_trials, p_contam, method, condition) |>
  print(n = Inf)

cat("\n=== Leading-edge depth moderation: sensitivity and drift error vs x_i ===\n\n")
bind_rows(
  results_long |>
    filter(process == process_labels["leading_edge"], p_contam > 0,
           condition != "All", method == "Flagging (ex-Gaussian)",
           measure %in% c("dep_cor_sens", "mcse_dep_cor_sens")),
  results_long |>
    filter(process == process_labels["leading_edge"], p_contam > 0,
           condition != "All",
           method %in% c("Simple", "Mixture (ex-Gaussian)",
                         "Threshold exclusion (p > .5)"),
           measure %in% c("dep_cor_edrift", "mcse_dep_cor_edrift"))
) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(n_trials, p_contam, method, condition, any_of(
    c("dep_cor_sens", "mcse_dep_cor_sens",
      "dep_cor_edrift", "mcse_dep_cor_edrift"))) |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  arrange(n_trials, p_contam, method, condition) |>
  print(n = Inf)

cat("\n=== Mixed contamination: pipelines (bias, contrast, correlations) ===\n\n")
results_long |>
  filter(process == process_labels["mixed"], p_contam > 0,
         method %in% pipeline_methods,
         measure %in% c("bias_drift", "bias_bound", "bias_ndt", "bias_effect",
                        "cor_drift", "cor_pcon")) |>
  pivot_wider(names_from = c(measure, condition), values_from = value) |>
  transmute(n_trials, p_contam, method,
            drift_easy = round(bias_drift_Easy, 3),
            drift_hard = round(bias_drift_Hard, 3),
            bound_easy = round(bias_bound_Easy, 3),
            ndt_easy   = round(bias_ndt_Easy, 3),
            effect     = round(bias_effect_All, 3),
            cor_drift_easy = round(cor_drift_Easy, 3),
            cor_pcon_easy  = round(cor_pcon_Easy, 3)) |>
  arrange(n_trials, p_contam, method) |>
  print(n = Inf)

cat("\n=== Mixed contamination: type-resolved flag sensitivity ===\n\n")
results_long |>
  filter(process == process_labels["mixed"], p_contam > 0,
         method == "Flagging (ex-Gaussian)",
         measure %in% c("clf_sens", "clf_sens_le", "clf_sens_dl",
                        "clf_sens_gu", "clf_spec", "clf_prec")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  select(n_trials, p_contam, condition, clf_sens, clf_sens_le, clf_sens_dl,
         clf_sens_gu, clf_spec, clf_prec) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  arrange(n_trials, p_contam, condition) |>
  print(n = Inf)
