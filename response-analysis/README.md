# Response analysis

This folder contains the corrected data, analysis scripts, results, and figures for the reply to Sánchez-Tójar and D'Amelio (2026) on Ruckman et al. (2024).

## Contents

- `data/` — Corrected dataset (`meta_complete_data2_corrected.csv`, 147 effect sizes from 72 studies). The corrections applied to the original dataset are:
  - `metafor::escalc(measure = "RBIS")` used for all biserial correlation calculations, fixing the parenthesis bug in the original custom rbis function.
  - Three missed effect sizes added (Carola et al. 2014 latency to first attack; Naretto and Chiaraviglio 2023 Rounds 2 and 3; Seaver and Hurd 2017 mirror-image aggression).
  - Two sample size corrections (Naretto and Chiaraviglio 2023 n = 48; Martin and Hengstebeck 1981 row 1 n = 19 from degrees of freedom).
  - Two sign corrections on reverse-coded measures (Zinzow-Kramer et al. 2015).
  - One effect-size substitution (Rose and Soole 2020, t from GLM in place of the omnibus F).
  - Two studies excluded (Yang et al. 2018; Podberscek and Serpell 1996) and eight rows dropped from Martin and Hengstebeck 1981 (rows 2 to 9).

- `scripts/` — Five R scripts that produce all results and figures reported in the reply.
  - `response_analysis_comparison.R` — main analysis. Runs six analyses on both the original and corrected datasets: (1) intercept-only meta-analysis, (2) publication year moderator, (3) effect-size source sensitivity, (4) publication bias via Egger's regression on residuals, (5) publication bias via Nakagawa et al. (2022) sqrt(1/n) moderator, (6) metafor::rma.mv sanity check with phylogenetic correlation matrix.
  - `no_fisher_z_check.R` — sensitivity analysis using raw rho in MCMCglmm (no Fisher Z transformation).
  - `metafor_raw_rho_check.R` — reproduces the exact framework used in Sánchez-Tójar and D'Amelio (2026) (raw r with metafor::rma.mv and phylogenetic correlation matrix).
  - `extend_all_methods.R` — extends the source-sensitivity and publication-bias analyses across all four combinations of transformation (raw r, Fisher Z) and framework (MCMCglmm, metafor::rma.mv).
  - `forest_plot_four_methods.R` — generates the four-method forest plot figure used in the reply.

- `results/` — All numeric outputs (CSV) referenced in the reply and its supplementary tables.

- `figures/` — Forest plot figures in PDF and PNG.

## Reproducibility

Set your R working directory to the root of this repository. Then run any script from `scripts/` directly. All paths inside the scripts are relative to the repository root.

For the full-length MCMC chains reported in the reply (Fisher Z + MCMCglmm and raw rho + MCMCglmm), each script uses 5,000,000 iterations with a 2,500,000 burn-in and a thinning interval of 1,000, yielding an effective sample size of approximately 2,500. Total runtime for the full chain-length version of `response_analysis_comparison.R` is 2 to 4 hours on a modern laptop. A `CHAIN_LENGTH <- "short"` toggle inside the script runs the same six analyses in about 15 minutes at the cost of noisier posteriors.

The phylogeny file (`Excel Sheets/list.nwk`) and the original dataset (`Excel Sheets/meta_complete_data2.csv`) are inherited from the main repository.

## Package versions

- `MCMCglmm` 2.36
- `metafor` 4.4 or later
- `meta` (for trimfill)
- `ape`, `phytools` (phylogeny handling)
- `ggplot2`, `patchwork` (figures)

## Release tag

This folder is tagged `response-analysis` on GitHub. The reply paper cites this exact release for reproducibility.
