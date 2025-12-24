# Metagenomic Spiked Primer Enrichment Analysis

R-based analysis scripts comparing MSSPE vs mNGS across three sequencing experiments. Each folder contains an R Markdown (`.Rmd`) file that renders publication-ready PDFs, figures, and tables.

## Structure
- `helper_utils.R`: Shared helpers (CSV processing, metadata parsing, export utilities, simulations).
- `data/`: Input data and per-experiment folders (`experiment_1/`, `experiment_2/`, `experiment_3/`, `experiment_all/`).
- Analysis markdown files (each has an `.Rmd` and outputs):
  	- `primer_tables/` (section 2.1)
	- `msspe_enrichment/` (section 2.2)
	- `sequencing_performance/` (section 2.2)
	- `sequencing_depth/` (section 2.3)
   	- `microbial_community/` (section 2.4)
- Outputs: `figures/`, `tables/`, and `*_files/` are auto-generated on knit or running all chunks.

## Requirements
- R ≥ 4.1 with XeLaTeX installed (TinyTeX or TeX Live).
- Packages: `tidyverse`, `dplyr`, `ggplot2`, `vegan`, `lme4`, `glmmTMB`, `patchwork`, `stringr`, `kableExtra`, `emmeans`, `taxize`, `magick`, `pdftools`, `superheat`.
