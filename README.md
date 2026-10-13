# Estimating Diffusion Model Parameters with the bmm R Package: A Tutorial

Companion repository for the tutorial paper on estimating Diffusion Decision Model (DDM) parameters using the [bmm](https://popov-lab.github.io/bmm/) R package. The paper demonstrates three DDM estimation approaches — the EZ-Diffusion Model, the Censored-Shifted Wald, and full DDM likelihood estimation — and shows that the hierarchical EZ-Diffusion Model evaluates condition effects on drift, boundary separation, and the relative starting point at a fraction of the cost of the trial-level approaches.

## Repository Structure

This GitHub repository contains the source code and data needed to reproduce all analyses. Rendered outputs (fitted models, figures, manuscript PDFs) are available on [OSF](https://osf.io/XXXXX/).

### Included in the GitHub repository

```text
├── manuscript/          Quarto manuscript source (.qmd) and references
│   ├── tutorial-ddm-bmm.qmd
│   └── references.bib
│
├── scripts/             Standalone R scripts, one per manuscript section
│   ├── 00_download_osf.R                     Download fitted models from OSF
│   ├── 01_prepare_lexical_data.R             Load, recode, and preprocess the lexical decision data
│   ├── 03_figure_diffusion_schematic.R       DDM schematic (Figure 1)
│   ├── 11_tutorial1_ezdm_3par.R              Part I, Tutorial 1: EZ-Diffusion 3-parameter version
│   ├── 12_tutorial2_ezdm_4par.R              Part I, Tutorial 2: EZ-Diffusion 4-parameter version
│   ├── 13_tutorial3_cswald.R                 Part II: Censored-Shifted Wald (simple + competing risks)
│   ├── 14_tutorial4_ddm.R                    Part II: full DDM likelihood estimation
│   ├── 15_cross_model_comparison.R           EZ-Diffusion frequency fits, Part III comparison figures
│   ├── 20_prepare_simulated_data.R           Generate the two simulated datasets
│   ├── 21_appendix_parameter_recovery.R      Appendix: recovery with an unbiased starting point
│   └── 22_bias_recovery.R                    Part I: recovery with a biased starting point
│
├── data/                Simulated and real datasets
│
└── functions/           Shared plotting and extraction utilities
    ├── clean_plot.R
    └── extract_fixef.R
```

### Available on OSF only

The following directories are excluded from the GitHub repository because they contain large binary files. They are archived on [OSF](https://osf.io/XXXXX/):

```text
├── output/              Cached model fits (.rds) — fitted brms/bmm objects
├── figures/             Generated figures (.pdf) — all manuscript figures
└── manuscript/*.pdf     Rendered manuscript (PDF, DOCX)
```

To reproduce locally, run the scripts in `scripts/` — fitted models will be cached in `output/` via the `file` argument in `bmm()`.

## Structure of the Tutorial

The manuscript centers on the hierarchical EZ-Diffusion Model as the primary analysis and is organized in three parts, preceded by a data-preprocessing section that introduces `bmm`'s contaminant-handling functions (`flag_contaminant_rts()`, `ezdm_summary_stats()`) on the lexical decision data of Wagenmakers et al. (2008).

**Part I — The EZ-Diffusion Model.** Tutorial 1 (`11_tutorial1_ezdm_3par.R`) fits the three-parameter version to the accuracy-coded data, tests condition effects on drift and bound, and refines the stimulus contrast to four word-frequency levels (`15_cross_model_comparison.R`, section 1). Tutorial 2 (`12_tutorial2_ezdm_4par.R`) switches to response coding and the four-parameter version to estimate the relative starting point. A recovery study (`22_bias_recovery.R`, data from `20_prepare_simulated_data.R`) shows what the four-parameter version recovers when the starting point is biased, against the three-parameter version and the full DDM.

**Part II — Trial-Level Estimation.** Tutorial 3 (`13_tutorial3_cswald.R`, `14_tutorial4_ddm.R`) fits the Censored-Shifted Wald (simple and competing-risks versions) and the full DDM at the same specification Part I arrived at.

**Part III — Comparison and Practical Guidelines.** `15_cross_model_comparison.R` compares the three approaches parameter by parameter and reports their computational cost; the manuscript turns the comparison into guidelines. It also refits the response-coded four-parameter EZ-Diffusion Model from mixture moments (to separate the aggregation-method offset from the approach) and applies the same per-cell 90%-interval check used for the EZ-Diffusion fits to the four trial-level fits, so fit quality is judged the same way for all three approaches. The appendix (`21_appendix_parameter_recovery.R`) reports parameter recovery with an unbiased starting point.

## Getting Started

Install required packages:

```r
install.packages("pacman")
pacman::p_load(here, bmm, brms, cmdstanr, tidyverse, tidybayes, patchwork, rtdists)
```

All analyses use `bmm` 1.4.0 from CRAN. The raw data are `speed_acc` from the `rtdists` package.

[cmdstanr](https://mc-stan.org/cmdstanr/) requires a working CmdStan installation. See `cmdstanr::install_cmdstan()`.

Open the RStudio project. To use pre-fitted models without re-running the analyses, first run `scripts/00_download_osf.R` to download cached model objects from OSF into `output/` (not yet functional: the OSF project is still to be set up). The tutorial scripts use the `file` argument in `bmm()`, so they will automatically load cached fits when available.

## References

- Ratcliff, R. (1978). A theory of memory retrieval. *Psychological Review, 85*(2), 59–108.
- Wagenmakers, E.-J., van der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. *Psychonomic Bulletin & Review, 14*(1), 3–22.
- Boag, R. J., et al. (2025). An expert guide to planning experimental tasks for evidence-accumulation modeling. *Advances in Methods and Practices in Psychological Science, 8*(2), 1–41.
- Frischkorn, G. T. & Popov, V. (2025). A tutorial for estimating Bayesian hierarchical mixture models for visual working memory tasks. *Behavior Research Methods*.
- bmm package: [https://popov-lab.github.io/bmm/](https://popov-lab.github.io/bmm/)

## OSF Repository

All data, code, and materials are also archived on OSF: [https://osf.io/XXXXX/](https://osf.io/XXXXX/)
