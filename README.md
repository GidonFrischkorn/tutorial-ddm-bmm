# Fitting Evidence Accumulation Models with the bmm R Package: A Tutorial

Companion repository for the tutorial paper on fitting evidence accumulation models (EAMs) using the [bmm](https://venpopov.github.io/bmm/) R package. The paper demonstrates three EAM implementations — EZ-Diffusion Model (EZDM), Censored-Shifted Wald Model (CSWALD), and Diffusion Decision Model (DDM) — with progressively complex tutorials.

## Repository Structure

This GitHub repository contains the source code and data needed to reproduce all analyses. Rendered outputs (fitted models, figures, manuscript PDFs) are available on [OSF](https://osf.io/XXXXX/).

### Included in the GitHub repository

```text
├── manuscript/          Quarto manuscript source (.qmd) and references
│   ├── tutorial-eam-bmm.qmd
│   └── references.bib
│
├── scripts/             Standalone R scripts for each tutorial
│   ├── prepare_simulated_data.R              Generate simulated datasets
│   ├── tutorial0_data_preprocessing.R        Data cleaning and summary stats
│   ├── tutorial1_ezdm_3par.R                 EZDM 3-parameter version
│   ├── tutorial2_ezdm_4par.R                 EZDM 4-parameter + model comparison
│   ├── tutorial3_cswald.R                    CSWALD simple + competing risks
│   ├── tutorial4_ddm.R                       Full DDM + cross-model comparison
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

**Tutorial 0 — The Impact of Contaminants** (`tutorial0_data_preprocessing.R`): Demonstrates empirically how contaminant RTs (fast guesses, slow outliers) bias EAM parameter estimates using simulated data with known ground truth. Introduces `flag_contaminant_rts()` for contaminant detection and `ezdm_summary_stats()` for computing summary statistics with simple, mixture, and robust methods. Verifies that preprocessing restores accurate parameter recovery.

**Tutorial 1 — EZDM (3-Parameter)** (`tutorial1_ezdm_3par.R`): Introduces the complete EAM workflow using the EZ-Diffusion Model with three parameters (drift, bound, ndt) on properly cleaned data. Covers model specification with `ezdm()`, fitting with `bmm()`, convergence diagnostics, posterior predictive checks, and hypothesis testing.

**Tutorial 2 — EZDM (4-Parameter)** (`tutorial2_ezdm_4par.R`): Extends the EZDM by adding the relative starting point (`zr`) parameter. Demonstrates how response biases can be estimated and compares the 3-parameter and 4-parameter versions using bridge sampling.

**Tutorial 3 — CSWALD as a Fast DDM Proxy** (`tutorial3_cswald.R`): Transitions to trial-level modeling with the Censored-Shifted Wald Model, positioned as a computationally efficient proxy for the DDM. Demonstrates both the simple version (censored approach for low-error tasks) and the competing-risks version. Compares parameter estimates against EZDM results from Tutorial 1.

**Tutorial 4 — DDM and Cross-Model Comparison** (`tutorial4_ddm.R`): Fits the full DDM using the Wiener first-passage time likelihood and systematically compares parameter recovery across all three model families (EZDM, CSWALD, DDM) against known ground truth. Includes a computational cost comparison, a real data application, and a practical decision tree for model selection.

## Getting Started

Install required packages:

```r
install.packages("pacman")
pacman::p_load(here, bmm, brms, cmdstanr, tidyverse, tidybayes, patchwork)
```

[cmdstanr](https://mc-stan.org/cmdstanr/) requires a working CmdStan installation. See `cmdstanr::install_cmdstan()`.

Open the RStudio project (`tutorial-eam-bmm.Rproj`). To use pre-fitted models without re-running the analyses, first run `scripts/00_download_osf.R` to download cached model objects from OSF into `output/`. The tutorial scripts use the `file` argument in `bmm()`, so they will automatically load cached fits when available.

## References

- Ratcliff, R. (1978). A theory of memory retrieval. *Psychological Review, 85*(2), 59–108.
- Wagenmakers, E.-J., van der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. *Psychonomic Bulletin & Review, 14*(1), 3–22.
- Boag, R. J., et al. (2025). An expert guide to planning experimental tasks for evidence-accumulation modeling. *Advances in Methods and Practices in Psychological Science, 8*(2), 1–41.
- Frischkorn, G. T. & Popov, V. (2025). A tutorial for estimating Bayesian hierarchical mixture models for visual working memory tasks. *Behavior Research Methods*.
- bmm package: [https://venpopov.github.io/bmm/](https://venpopov.github.io/bmm/)

## OSF Repository

All data, code, and materials are also archived on OSF: [https://osf.io/XXXXX/](https://osf.io/XXXXX/)
