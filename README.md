# metagenomic-spiked-primer-enrichment-paper
scripts used in the publication

## Export helper functions

Analysis notebooks can now call shared helpers from `helper_utils.R` to save publication-ready assets without re-implementing boilerplate:

- `ensure_export_dirs(fig_dir, table_dir)` creates the figure/table directories if they do not exist.
- `export_plot_png(plot_obj, fig_dir, filename, width, height, dpi)` writes high-resolution PNGs from any ggplot object.
- `export_table_csv(table_df, table_dir, filename)` saves companion CSV tables via `readr::write_csv()`.
- `export_table_png(latex_table, table_dir, output_basename, paper_size, density)` compiles a LaTeX table (e.g., a `kableExtra` object) to PDF with TinyTeX and converts it to a 600 dpi PNG using `magick`/`pdftools`.

See `msspe_enrichment/MSSPE_enrichment.Rmd` for an example of wiring these helpers into an export chunk; other R Markdown files can reuse the same pattern by sourcing `helper_utils.R`.
