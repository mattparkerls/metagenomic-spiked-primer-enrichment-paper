# Metagenomic Spiked Primer Enrichment Analysis

R-based analysis modules comparing MSSPE vs NGS across three sequencing experiments. Each module is an R Markdown (`.Rmd`) that renders publication-ready PDFs, figures, and tables.

## Structure
- `helper_utils.R`: Shared helpers (CSV processing, metadata parsing, export utilities, simulations).
- `data/`: Input data and per-experiment folders (`experiment_1/`, `experiment_2/`, `experiment_3/`, `experiment_all/`).
- Analysis modules (each has an `.Rmd` and outputs):
	- `msspe_enrichment/`
	- `sequencing_performance/`
	- `microbial_community/`
	- `sequencing_depth/`
	- `pathogen_investigation/`
	- `primer_detection/`
	- `msspe_regression/`
- Outputs: `figures/`, `tables/`, and `*_files/` are auto-generated on knit.

## Requirements
- R ≥ 4.1 with XeLaTeX installed (TinyTeX or TeX Live).
- Packages: `tidyverse`, `dplyr`, `ggplot2`, `vegan`, `lme4`, `glmmTMB`, `patchwork`, `stringr`, `kableExtra`, `emmeans`, `taxize`, `magick`, `pdftools`.
