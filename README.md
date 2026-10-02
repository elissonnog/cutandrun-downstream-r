# CUT&RUN downstream analysis in R

A compact, parameterized workflow for downstream analysis of CUT&RUN consensus-peak counts. It accepts an explicit peak-by-sample integer count matrix, sample metadata, an additive categorical design, and declared contrasts. It produces sample QC, a DESeq2 differential-peak analysis, normalized counts, figures, result tables, and an HTML report.

This workflow model was developed during my postdoctoral research at Van Andel Institute. This reusable packaging is distinct from the public study reanalysis below and from the original VAI study code. The public example demonstrates reproducible mechanics; it does not claim authorship of the source data or its biological conclusions.

## Public GSE186608 reanalysis

The primary example reanalyzes six HCT116 H3K27ac CUT&RUN samples from [Ivancevic et al., *Science Advances* (2024), DOI 10.1126/sciadv.ado1218](https://doi.org/10.1126/sciadv.ado1218): untreated, cobimetinib, and TNFalpha, with two treatment replicates per condition. Source files are pinned to tag `v1.0`, commit `faeef682bebdfb982c70b47a5d7df4ac77f61264`, and verified by SHA-256 before use.

![Cobimetinib versus untreated MA plot](figures/gse186608-cobimetinib-ma.png)

*GSE186608 public-data reanalysis. MA plot of unshrunken DESeq2 estimates for cobimetinib versus untreated; color denotes adjusted-p-value significance at FDR 0.05.*

![Cobimetinib significant-peak replicate heatmap](figures/gse186608-cobimetinib-significant-heatmap.png)

*GSE186608 public-data reanalysis. Row-centered VST values for the deterministically ranked top 30 of 2,743 significant cobimetinib-associated regions. Samples remain in input order; the top bar and legend show condition.*

All 38,649 input regions remain in exported result tables; 31,730 nonzero regions receive a DESeq2 test statistic. At FDR 0.05, cobimetinib versus untreated has 2,743 significant regions (1,084 increased and 1,659 decreased relative signals), while TNFalpha versus untreated has 914 (802 increased and 112 decreased).

These are limited claims. The model is `~ condition`, with no invented batch or pairing term, and only two treatment replicates per condition. Untreated replicate 2 was resequenced and merged by the source authors and remains one sample. Counts are `bedtools multicov` aligned-read overlaps, not paired-fragment counts, and the region set is the union of per-sample top-20k peaks. No IgG/background or spike-in normalization was reported for these matrices, so DESeq2 estimates relative signal only; it does not establish absolute or global H3K27ac gain or loss. Raw-read and BAM-level QC are outside this repository.

## Requirements

- R
- Pandoc
- CRAN: `rmarkdown`, `knitr`, `ggplot2`, `openssl`
- Bioconductor: `DESeq2`, `SummarizedExperiment`

Optional coordinate annotation is a functioning, guarded stage rather than a placeholder. For the source-supported human hg38 path it requires these official Bioconductor packages:

- `GenomicRanges`
- `IRanges`
- `GenomeInfoDb`
- `AnnotationDbi`
- `ChIPseeker`
- `TxDb.Hsapiens.UCSC.hg38.knownGene`
- `org.Hs.eg.db`

Annotation additionally requires real hg38 peak coordinates with UCSC-style sequence names compatible with the selected TxDb.

The reviewed local runtime used R 4.6.1, Bioconductor 3.23, Pandoc 3.12, `DESeq2` 1.52.0, `GenomicRanges` 1.64.0, `IRanges` 2.46.0, `GenomeInfoDb` 1.48.0, `AnnotationDbi` 1.74.0, `ChIPseeker` 1.48.0, `TxDb.Hsapiens.UCSC.hg38.knownGene` 3.22.0, and `org.Hs.eg.db` 3.23.1.

The workflow never installs packages at run time.

## Inputs

- `peak_counts.csv`: required columns `peak_id`, `seqnames`, `start`, and `end`, followed by non-negative integer sample columns. Peak IDs and coordinate intervals must be unique and coordinates must already describe a consensus peak set.
- `metadata.csv`: one row per sample, with `sample_id` plus categorical design variables.
- `contrasts.csv`: columns `label`, `factor`, `numerator`, and `denominator`.
- An R configuration file modeled on `config/gse186608_config.R`, including explicit `genome_build` and `coordinate_system` fields.

Supported coordinate declarations are `one_based_closed` and `bed_zero_based_half_open`. BED starts are converted once to one-based closed coordinates (`start + 1`, unchanged end); all result and annotation tables use one-based closed coordinates. Other conventions fail rather than being guessed.

Validation covers blank or duplicate headers (case-insensitive), unique IDs and intervals, integer/nonnegative count semantics, declared coordinate validity and conversion, exact sample matching and safe reordering, additive categorical designs, full-rank model matrices, declared contrast levels, unequal numerator/denominator levels, and case-insensitive output-name collisions.

## Genomic annotation

When annotation is enabled, the workflow:

1. requires the annotation and peak `genome_build` declarations to match exactly;
2. loads the configured TxDb and verifies its own genome metadata against that build;
3. requires every peak seqname to exist in the TxDb;
4. builds `GRanges` from the internally one-based closed peak coordinates;
5. runs `ChIPseeker::annotatePeak` with the configured TxDb, OrgDb, and TSS window;
6. preserves and verifies every unique `peak_id`; and
7. adds annotation columns to each differential-peak table as well as `peak_annotation.csv`.

The implemented annotation path supports human hg38 and refuses genome-build mixing. Peak coordinates, the configured build, and the selected TxDb must agree before annotation runs. This resolves an inconsistency in the archived VAI scripts, where hg38 variables coexisted with a separate hg19 gene-overlap block.

The core synthetic fixture uses `genome_build = "synthetic_v1"` and `chrSynthetic`; annotation is intentionally disabled. Simply toggling annotation on fails because the configured hg38 annotation build does not match the synthetic peak build. The separate `synthetic_hg38_annotation` mechanics fixture uses public hg38 TxDb-derived intervals with synthetic counts; it contains no experimental observations and supports annotation-path testing without relabeling synthetic coordinates.

## Normalization choice

The default `deseq2_median_ratio` mode estimates relative abundance and assumes that most retained peaks do not shift together. This can be inappropriate when CUT&RUN has a genuine global occupancy change. Such a run must not be interpreted as proving global gain or loss.

If an upstream experimental strategy provides defensible spike-in or other external scaling factors, choose `provided_deseq2_size_factors` and supply a two-column `sample_id,deseq2_size_factor` CSV. The factor is explicitly a DESeq2 divisor: `normalized count = count / deseq2_size_factor`. A multiplicative spike-in scaling factor must be inverted before use and its provenance documented. The workflow validates positivity, sample identity, reordering, geometric-mean centering, and final assignment, but it cannot establish whether the upstream control or transformation is biologically valid.

## Reproduce the public example

From the repository root:

```sh
Rscript --vanilla scripts/prepare_gse186608.R
Rscript --vanilla tests/test_helpers.R
Rscript --vanilla tests/test_core_analysis.R
Rscript --vanilla tests/test_annotation.R
Rscript --vanilla run_analysis.R --check-dependencies
Rscript --vanilla run_analysis.R config/gse186608_config.R
```

Review the saved [standalone HTML report](report/gse186608_h3k27ac_downstream.html), [contrast summary](results/gse186608/contrast_summary.csv), and [top significant-region snapshot](results/gse186608/top_significant_regions.csv). Complete generated tables remain under the ignored `output/gse186608_h3k27ac/` directory after a local run.

The public files come from the authors' [Zenodo record 10996183](https://zenodo.org/records/10996183), licensed CC BY 4.0. This workflow verifies the pinned files, joins the exactly matching BED coordinates to count rows, adds explicit GEO sample metadata and contrasts, converts BED coordinates for Bioconductor, runs DESeq2, and adds hg38 annotation and figures. Source files are not committed; pins and checksums are in `config/gse186608_sources.csv`. This attribution does not grant a new license for the workflow code.

Synthetic negative-binomial fixtures remain secondary regression tests. Regenerate them with:

```sh
Rscript --vanilla example/synthetic/generate_fixture.R
Rscript --vanilla example/synthetic_hg38_annotation/generate_fixture.R
```

## Outputs

The configured output directory receives:

- `cutandrun_downstream.html`
- `tables/run_metadata.csv`
- `tables/normalized_peak_counts.csv`
- `tables/size_factors.csv`
- `tables/contrast_summary.csv`
- `tables/<contrast>_deseq2.csv`
- `tables/peak_annotation.csv` when annotation succeeds; the same annotation columns are joined to every contrast result
- `figures/library_qc.png`
- `figures/pca.png`
- `figures/sample_correlation.png`
- `figures/<contrast>_ma.png`
- `figures/<contrast>_volcano.png`
- `figures/<contrast>_significant_peak_heatmap.png`

## Scope and provenance

The workflow begins with a consensus peak-by-sample count matrix. Upstream alignment, BAM processing, consensus-peak definition, spike-in factor derivation, and raw-read QC remain outside this repository.

The reusable implementation focuses on differential-peak summaries and transformed-count QC and heatmaps. It does not claim to implement DiffBind, time-course clustering, cross-mark integration, regulatory networks, or pathway enrichment. External size factors are accepted only when a defensible upstream experimental strategy has already produced them; the workflow validates and applies those factors but does not infer their biological provenance.

No repository license has been selected or added.
