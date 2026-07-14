# Assessing the association between animal color and behavior

Data and analysis code for a phylogenetic Bayesian meta-analysis of the association between color and aggression across animal taxa, and for a subsequent reply to a published critique.

## Papers

**Original paper.** Ruckman, S. N., Humphrey, E. A., Muzzey, L., Prantalou, I., Pleasants, M., and Hughes, K. A. (2024). Assessing the association between animal color and behavior: a meta-analysis of experimental studies. *Ecology and Evolution* 14: e70655.

**Reply paper.** Ruckman, S. N. and Humphrey, E. A. Toward shared standards in evolutionary meta-analysis: a reply to Sánchez-Tójar and D'Amelio (2026). *Ecology and Evolution* (in press). DOI: [to be added]

## Repository structure

### `Excel Sheets/`

Data files used in the original meta-analysis and its supporting searches.

- `meta_complete_data2.csv` — original compiled dataset (169 effect sizes from 74 studies).
- `list.nwk` — phylogenetic tree in Newick format.
- `species_order.csv`, `species list.txt` — species metadata.
- `scopus search.xlsx`, `wos updated search.xlsx` — records of the literature searches used to identify included studies.
- `Papers_with_color_and_aggression.xlsx` — full list of screened papers with inclusion decisions.
- `obsolete/` — earlier versions of the data files and scripts, retained for provenance.

### `R Files/`

Analysis code and saved models for the original paper.

- `Effect Size Calculations.R` — code to compute effect sizes from the raw study data.
- `MCMCglmm_Z.R` — main MCMCglmm meta-analysis on Fisher-Z transformed correlations.
- `agg_Z.RDATA`, `allRnd_Z.RDATA`, `class_Z.RDATA`, `rsocial_Z.RDATA`, `rsocolor_Z.RDATA` — saved MCMCglmm model objects.
- `Figures/` — final figures used in the manuscript, and supporting versions.
- `obsolete/` — earlier scripts and model objects, retained for provenance.

### `response-analysis/`

Corrected data, analysis code, results, and figures for the reply paper. This folder is self-contained. See `response-analysis/README.md` for details.

Tagged as the release `response-analysis` on GitHub.

## Reproducibility

Set your R working directory to the root of this repository. All scripts in `R Files/` and `response-analysis/scripts/` use paths relative to the repo root.

### Package versions

- `MCMCglmm` 2.36
- `metafor` 4.4 or later
- `meta` (for trim-and-fill)
- `ape`, `phytools` (phylogeny handling)
- `ggplot2`, `patchwork` (figures)

### Runtime

The full-length MCMC chains used for both the original paper and the reply run for 5,000,000 iterations with a 2,500,000 burn-in and a thinning interval of 1,000. Total runtime for the main comparison script (`response-analysis/scripts/response_analysis_comparison.R`) is 2 to 4 hours on a modern laptop. A `CHAIN_LENGTH <- "short"` toggle inside the script runs the same six analyses in about 15 minutes.

## Citation

If you use data or code from this repository, please cite the original paper (Ruckman et al. 2024). If you use the corrected dataset or any analysis in `response-analysis/`, please also cite the reply paper.

## Contact

Sarah N. Ruckman, University of California, Irvine. sarahnruckman@gmail.com
