# Estimating Diffusion Model Parameters with the bmm R Package: A Tutorial

Companion repository for the tutorial paper on estimating Diffusion Decision Model (DDM) parameters with the [bmm](https://popov-lab.github.io/bmm/) R package. The paper demonstrates three DDM estimation approaches, the EZ-Diffusion Model, the Censored-Shifted Wald, and full DDM likelihood estimation, on the lexical decision data of Wagenmakers et al. (2008), and makes the case for the hierarchical EZ-Diffusion Model as a primary analysis of condition effects on drift, boundary separation, non-decision time, and the relative starting point.

- Manuscript, Supplementary Material, fitted models, and figures: [OSF](https://osf.io/7y54v/)
- Source code and data: this repository

## Repository Structure

```text
├── manuscript/
│   ├── tutorial-ddm-bmm.qmd       Quarto manuscript (apaquarto format)
│   ├── supplement.qmd             Supplementary Material S1: parameter recovery
│   ├── references.bib
│   └── _extensions/               apaquarto extension (pinned)
│
├── scripts/                       Standalone R scripts, one per manuscript section
│   ├── 00_download_osf.R                 Download fitted models and figures from OSF
│   ├── 01_prepare_lexical_data.R         Load, recode, and preprocess the lexical decision data
│   ├── 03_figure_diffusion_schematic.R   DDM schematic (Figure 1)
│   ├── 11_tutorial1_ezdm_3par.R          Part I, Tutorial 1: EZ-Diffusion 3-parameter version
│   ├── 12_tutorial2_ezdm_4par.R          Part I, Tutorial 2: EZ-Diffusion 4-parameter version
│   ├── 13_tutorial3_cswald.R             Part II: Censored-Shifted Wald (simple + competing risks)
│   ├── 14_tutorial3_ddm.R                Part II: full DDM likelihood estimation
│   ├── 15_cross_model_comparison.R       EZ-Diffusion frequency fits, Part III comparison figures
│   ├── 20_prepare_simulated_data.R       Generate the two simulated datasets
│   ├── 21_appendix_parameter_recovery.R  Supplement S1: recovery with an unbiased starting point
│   ├── 22_bias_recovery.R                Part I: recovery with a biased starting point
│   └── 99_upload_osf.R                   Maintainer: upload fits, figures, and PDFs to OSF
│
├── data/                          Prepared lexical decision data and simulated datasets (CSV)
├── functions/                     Shared plotting and extraction utilities
├── output/                        Fitted model objects (.rds), from OSF or from the scripts
└── figures/                       Manuscript figures (.pdf), from OSF or from the scripts
```

`output/` and `figures/` are not in the repository because of their size (about 250 MB of fitted models). `scripts/00_download_osf.R` downloads them from OSF; alternatively, the scripts regenerate them.

## Structure of the Tutorial

The manuscript is organized in three parts, preceded by a data-preprocessing section that introduces `bmm`'s contaminant-handling functions (`flag_contaminant_rts()`, `ezdm_summary_stats()`).

**Part I: The EZ-Diffusion Model.** Tutorial 1 (`11_tutorial1_ezdm_3par.R`) fits the three-parameter version to the accuracy-coded data, tests condition effects on drift and bound, and refines the stimulus contrast to four word-frequency levels (`15_cross_model_comparison.R`, section 1). Tutorial 2 (`12_tutorial2_ezdm_4par.R`) switches to response coding and the four-parameter version to estimate the relative starting point. A recovery study (`22_bias_recovery.R`, data from `20_prepare_simulated_data.R`) shows what the four-parameter version recovers when the starting point is biased.

**Part II: Trial-Level Estimation.** Tutorial 3 (`13_tutorial3_cswald.R`, `14_tutorial3_ddm.R`) fits the Censored-Shifted Wald (simple and competing-risks versions) and the full DDM at the specification Part I arrived at.

**Part III: Comparison and Practical Guidelines.** `15_cross_model_comparison.R` compares the three approaches parameter by parameter, reports their computational cost, and applies the same per-cell posterior predictive check to all of them. The Supplementary Material (`21_appendix_parameter_recovery.R`) reports parameter recovery with an unbiased starting point.

## Reproducing the Analyses

Software versions used for the reported fits and the rendered manuscript:

| Software | Version |
|---|---|
| R | 4.6.1 |
| bmm | 1.4.1 (CRAN) |
| brms | 2.23.0 |
| cmdstanr | 0.9.0 |
| CmdStan | 2.40.0 |

Install the packages (cmdstanr needs a working CmdStan; see `cmdstanr::install_cmdstan()`):

```r
install.packages("pacman")
pacman::p_load(here, bmm, brms, cmdstanr, tidyverse, tidybayes, patchwork, rtdists, osfr)
```

Then open the RStudio project and either

1. run `scripts/00_download_osf.R` to download the fitted models and figures (the scripts load cached fits through the `file` argument of `bmm()` and skip sampling), or
2. run the scripts in this order to refit everything: `01` → `20` → `11` → `12` → `13` → `14` → `15` → `21` → `22`.

Fitting times on an Apple M3 Max (five chains, two threads per chain, one fit at a time): the EZ-Diffusion fits take seconds to minutes; the trial-level fits on the lexical decision data take 3 to 12 hours each; the DDM fits on the simulated data take 1 to 4 hours. All fits use `backend = "cmdstanr"`.

The manuscript is rendered with Quarto (`quarto render manuscript/tutorial-ddm-bmm.qmd`); every number in the text is computed from the fitted models in `output/` and the CSVs in `data/` when the document is rendered.

## Data

The lexical decision data are `speed_acc` from the `rtdists` package (Singmann et al., 2024, version 0.11-5; GPL-3), originally reported by Wagenmakers, Ratcliff, Gomez, and McKoon (2008). The CSVs in `data/` are the prepared versions produced by `01_prepare_lexical_data.R` (contaminant probabilities and summary statistics added) and the simulated datasets produced by `20_prepare_simulated_data.R`.

## License

The code in this repository is released under the MIT License (see `LICENSE`). The lexical decision data remain subject to the terms of the `rtdists` package.

## References

- Ratcliff, R. (1978). A theory of memory retrieval. *Psychological Review, 85*(2), 59–108.
- Wagenmakers, E.-J., van der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. *Psychonomic Bulletin & Review, 14*(1), 3–22.
- Wagenmakers, E.-J., Ratcliff, R., Gomez, P., & McKoon, G. (2008). A diffusion model account of criterion shifts in the lexical decision task. *Journal of Memory and Language, 58*(1), 140–159.
- Boag, R. J., et al. (2025). An expert guide to planning experimental tasks for evidence-accumulation modeling. *Advances in Methods and Practices in Psychological Science, 8*(2), 1–41.
- Frischkorn, G. T., & Popov, V. (2025). A tutorial for estimating Bayesian hierarchical mixture models for visual working memory tasks. *Behavior Research Methods, 57*. https://doi.org/10.3758/s13428-025-02643-0
- Singmann, H., Brown, S., Gretton, M., & Heathcote, A. (2024). rtdists: Response time distributions (R package version 0.11-5). https://CRAN.R-project.org/package=rtdists
- bmm package: [https://popov-lab.github.io/bmm/](https://popov-lab.github.io/bmm/)
