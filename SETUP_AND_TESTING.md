# Setup, custom data, and testing

Run every command from the repository root. Configuration paths may be absolute; relative paths are resolved from the repository root.

## Install requirements

Install a current R release and external Pandoc. Confirm both executables are on `PATH`:

```sh
R --version
pandoc --version
```

Pandoc installation instructions and binaries are available from the [official Pandoc site](https://pandoc.org/installing.html). Common package-manager commands are `brew install pandoc` on macOS and `sudo apt install pandoc` on Debian/Ubuntu.

Install the core R packages:

```r
install.packages(c("BiocManager", "rmarkdown", "knitr", "ggplot2", "openssl"))
BiocManager::install(c("DESeq2", "SummarizedExperiment"))
```

Genomic annotation is optional. When enabled, install the annotation framework and the exact TxDb/OrgDb packages declared in the configuration. For the included hg38 example:

```r
BiocManager::install(c(
  "GenomicRanges", "IRanges", "GenomeInfoDb", "AnnotationDbi", "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"
))
```

No workflow command installs packages automatically. The last verified environment used R 4.6.1, Pandoc 3.12, rmarkdown 2.32, knitr 1.52, ggplot2 4.0.3, openssl 2.4.2, DESeq2 1.52.0, SummarizedExperiment 1.42.0, GenomicRanges 1.64.0, IRanges 2.46.0, GenomeInfoDb 1.48.0, AnnotationDbi 1.74.0, ChIPseeker 1.48.0, TxDb.Hsapiens.UCSC.hg38.knownGene 3.22.0, and org.Hs.eg.db 3.23.1. This is a tested-version record, not a lockfile.

## Custom-data contract

Use three CSV files:

- Peak counts: columns must be `peak_id,seqnames,start,end` followed only by sample columns. Peak IDs, intervals, and headers must be unique and nonblank. Coordinates and counts must be finite integers; counts must be nonmissing and nonnegative, with at least two samples. Unexpected annotation or other columns are rejected rather than treated as counts.
- Metadata: must contain a unique, nonblank `sample_id` column. Its values must exactly match the count-column names; order may differ and is realigned safely. Design variables must be complete categorical values with at least two levels and descriptive, nonnumeric labels.
- Contrasts: must contain `label,factor,numerator,denominator`. Labels must be unique and not collide after filename sanitization. `factor` must occur in the design, numerator and denominator must differ, and both levels must exist in metadata.

Only additive categorical designs such as `~ condition` or `~ batch + condition` are supported. Interactions, transformations, continuous covariates, response variables, and rank-deficient designs are rejected.

Choose one coordinate convention and declare it exactly:

- `one_based_closed`: requires `start >= 1` and `end >= start`; coordinates are used as supplied.
- `bed_zero_based_half_open`: requires `start >= 0` and `end > start`; start is converted once to one-based closed coordinates for analysis.

The synthetic files are small schema examples:

- [peak_counts.csv](example/synthetic/peak_counts.csv)
- [metadata.csv](example/synthetic/metadata.csv)
- [contrasts.csv](example/synthetic/contrasts.csv)

## Configure and run custom data

1. Copy the complete example configuration:

```sh
cp config/example_config.R config/my_analysis.R
```

2. Edit `config/my_analysis.R`:

   - Set `project_title`, `data_label`, the three input paths, and a new `output_dir`.
   - Declare `coordinate_system`, the exact `genome_build`, additive `design`, `min_total_count`, `fdr`, and plotting limits.
   - Choose `normalization$method = "deseq2_median_ratio"` for relative signal, or `"provided_deseq2_size_factors"` with `size_factors_path`.
   - Keep annotation disabled for synthetic or unknown coordinates. To enable it, set `annotation$enabled = TRUE`; `annotation$genome_build` must exactly equal `genome_build`, the selected TxDb must declare that build and compatible seqnames, and `tss_region` must span zero.

For provided size factors, the CSV must contain exactly `sample_id,deseq2_size_factor`, cover the same samples once each, and contain finite positive values. Values are geometrically centered and used as DESeq2 divisors. The workflow validates the file but cannot validate the biological derivation or spike-in convention.

3. Check the exact configuration. A missing required dependency produces a nonzero exit and nothing is installed:

```sh
Rscript --vanilla run_analysis.R --check-dependencies config/my_analysis.R
```

4. Render the analysis:

```sh
Rscript --vanilla run_analysis.R config/my_analysis.R
```

The configured output directory receives the HTML report, full result ledgers, normalized peak counts, size factors, diagnostics, optional annotation and coverage, input/code hashes, resolved configuration, software versions, warnings, and session information.

## Tests and public-example regeneration

```sh
Rscript --vanilla tests/test_helpers.R
Rscript --vanilla tests/test_core_analysis.R
Rscript --vanilla tests/test_annotation.R
Rscript --vanilla tests/test_preflight.R
```

To regenerate the public report and compact portfolio files:

```sh
Rscript --vanilla run_analysis.R config/gse186608_config.R
Rscript --vanilla scripts/export_portfolio_results.R
Rscript --vanilla scripts/render_portfolio_snapshot.R
```
