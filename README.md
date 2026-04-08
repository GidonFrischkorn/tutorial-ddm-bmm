# Estimating Diffusion Model Parameters with the bmm R Package: A Tutorial

Companion repository for the tutorial paper on estimating Diffusion Decision Model (DDM) parameters using the [bmm](https://venpopov.github.io/bmm/) R package. The paper demonstrates three DDM estimation approaches — the EZ-Diffusion Model, the Censored-Shifted Wald, and full DDM likelihood estimation — with progressively complex tutorials.

## Repository Structure

This GitHub repository contains the source code and data needed to reproduce all analyses. Rendered outputs (fitted models, figures, manuscript PDFs) are available on [OSF](https://osf.io/XXXXX/).

### Included in the GitHub repository

```text
├── manuscript/          Quarto manuscript source (.qmd) and references
│   ├── tutorial-ddm-bmm.qmd
│   └── references.bib
│
├── scripts/             Standalone R scripts for each tutorial
│   ├── prepare_simulated_data.R              Generate simulated datasets
│   ├── tutorial0_data_preprocessing.R        Data cleaning and summary stats
│   ├── tutorial1_ezdm_3par.R                 EZ-Diffusion 3-parameter version
│   ├── tutorial2_ezdm_4par.R                 EZ-Diffusion 4-parameter + model comparison
│   ├── tutorial3_cswald.R                    Censored-Shifted Wald simple + competing risks
│   ├── tutorial4_ddm.R                       Full DDM + cross-approach comparison
│   └── 00_download_osf.R                     Download fitted models from OSF
│
├── data/                Simulated and real datasets
│
└── functions/           Shared plotting utilities
    └── clean_plot.R
```

### Available on OSF only

The following directories are excluded from the GitHub repository because they contain large binary files. They are archived on [OSF](https://osf.io/XXXXX/):

```text
├── output/              Cached model fits (.rds) — fitted brms/bmm objects
├── figures/             Generated figures (.pdf) — all manuscript figures
└── manuscript/*.pdf     Rendered manuscript (PDF, DOCX)
```

To reproduce locally, run the scripts in `scripts/` — fitted models will be cached in `output/` via the `file` argument in `bmm()`.

## Tutorials

**Tutorial 0 — The Impact of Contaminants** (`tutorial0_data_preprocessing.R`): Demonstrates empirically how contaminant RTs (fast guesses, slow outliers) bias DDM parameter estimates using simulated data with known ground truth. Introduces `flag_contaminant_rts()` for contaminant detection and `ezdm_summary_stats()` for computing summary statistics with simple, mixture, and robust methods. Verifies that preprocessing restores accurate parameter recovery.

**Tutorial 1 — EZ-Diffusion Model (3-Parameter)** (`tutorial1_ezdm_3par.R`): Introduces the complete DDM estimation workflow using the EZ-Diffusion Model with three parameters (drift, bound, ndt) on properly cleaned data. Covers model specification with `ezdm()`, fitting with `bmm()`, convergence diagnostics, posterior predictive checks, and hypothesis testing.

**Tutorial 2 — EZ-Diffusion Model (4-Parameter)** (`tutorial2_ezdm_4par.R`): Extends the EZ-Diffusion Model by adding the relative starting point (`zr`) parameter. Demonstrates how response biases can be estimated and compares the 3-parameter and 4-parameter versions.

**Tutorial 3 — Censored-Shifted Wald** (`tutorial3_cswald.R`): Transitions to trial-level modeling with the Censored-Shifted Wald, a computationally efficient trial-level approach for estimating DDM parameters. Demonstrates both the simple version (censored approach for low-error tasks) and the competing-risks version. Compares parameter estimates against the EZ-Diffusion results from Tutorial 1.

**Tutorial 4 — Full DDM Likelihood Estimation** (`tutorial4_ddm.R`): Fits the full DDM using the Wiener first-passage time likelihood and systematically compares parameter recovery across all three estimation approaches. Includes a computational cost comparison, a real data application, and a practical decision tree for approach selection.

## Getting Started

Install required packages:

```r
install.packages("pacman")
pacman::p_load(here, bmm, brms, cmdstanr, tidyverse, tidybayes, patchwork)
```

[cmdstanr](https://mc-stan.org/cmdstanr/) requires a working CmdStan installation. See `cmdstanr::install_cmdstan()`.

Open the RStudio project. To use pre-fitted models without re-running the analyses, first run `scripts/00_download_osf.R` to download cached model objects from OSF into `output/`. The tutorial scripts use the `file` argument in `bmm()`, so they will automatically load cached fits when available.

## References

- Ratcliff, R. (1978). A theory of memory retrieval. *Psychological Review, 85*(2), 59–108.
- Wagenmakers, E.-J., van der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. *Psychonomic Bulletin & Review, 14*(1), 3–22.
- Boag, R. J., et al. (2025). An expert guide to planning experimental tasks for evidence-accumulation modeling. *Advances in Methods and Practices in Psychological Science, 8*(2), 1–41.
- Frischkorn, G. T. & Popov, V. (2025). A tutorial for estimating Bayesian hierarchical mixture models for visual working memory tasks. *Behavior Research Methods*.
- bmm package: [https://venpopov.github.io/bmm/](https://venpopov.github.io/bmm/)

## OSF Repository

All data, code, and materials are also archived on OSF: [https://osf.io/XXXXX/](https://osf.io/XXXXX/)
