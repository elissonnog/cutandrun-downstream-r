# CUT&RUN downstream analysis in R

An R workflow for quality control, DESeq2 differential peak analysis, hg38 annotation, replicate-level visualization, and a reviewable HTML report from a consensus peak-count matrix. The workflow pattern was developed during my postdoctoral work at Van Andel Institute and is demonstrated here with public HCT116 H3K27ac data.

## Dataset and design

The main example uses H3K27ac CUT&RUN data from [GSE186608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186608), reported by Ivancevic et al. (*Science Advances*, 2024). It contains HCT116 cells under three 24-hour conditions, with two biological replicates each: untreated, 1 uM cobimetinib, and 100 ng/mL TNFalpha.

- Comparisons: cobimetinib versus untreated; TNFalpha versus untreated
- Model: `~ condition`
- Input: aligned-read overlaps across the authors' union of per-sample top-20k peaks
- Genome/coordinates: hg38; source BED converted once from zero-based half-open to one-based closed
- Significance: Benjamini-Hochberg adjusted p-value < 0.05

The design follows the available sample metadata. The authors' resequenced/merged untreated replicate is represented as one biological sample.

## Run the public example

From the repository root:

```sh
Rscript --vanilla scripts/prepare_gse186608.R
Rscript --vanilla run_analysis.R config/gse186608_config.R
```

The preparation script downloads the two author-supplied source files, verifies their SHA-256 checksums and dimensions, and validates the exact BED/count-row mapping. The analysis writes complete tables and figures under `output/gse186608_h3k27ac/`.

Review the saved [HTML report](report/gse186608_h3k27ac_downstream.html), [contrast summary](results/gse186608/contrast_summary.csv), and [top significant-region snapshot](results/gse186608/top_significant_regions.csv). The workflow requires R, Pandoc, `rmarkdown`, `knitr`, `ggplot2`, `openssl`, `DESeq2`, `SummarizedExperiment`, and the hg38 Bioconductor annotation packages named in the configuration.

Dependency checks, tested package versions, tests, and result-regeneration commands are documented in [SETUP_AND_TESTING.md](SETUP_AND_TESTING.md).

## Results

All 38,649 input regions are retained in exported result ledgers. Of these, 6,919 are all-zero and 31,730 receive a finite raw p-value; DESeq2 independent filtering leaves 26,257 finite adjusted p-values for cobimetinib and 24,433 for TNFalpha. At FDR 0.05:

- Cobimetinib versus untreated: 2,743 significant regions (1,084 increased; 1,659 decreased).
- TNFalpha versus untreated: 914 significant regions (802 increased; 112 decreased).

![Cobimetinib versus untreated MA plot](figures/gse186608-cobimetinib-ma.png)

*Mean H3K27ac count versus unshrunken DESeq2 log2 fold change; colored points pass FDR 0.05.*

![Heatmap of cobimetinib-associated regions](figures/gse186608-cobimetinib-significant-heatmap.png)

*Row-centered variance-stabilized counts for the 30 highest-ranked significant regions, with both replicates per condition shown.*

## Comparison with the published study

Ivancevic et al. report that MAPK/AP1 inhibition by cobimetinib preferentially decreases H3K27ac at responsive enhancers, whereas TNFalpha stimulation increases it, including at LTR10-derived enhancers. The aggregate significant-region counts are qualitatively consistent with these treatment-specific results: losses outnumber gains after cobimetinib (1,659 versus 1,084), whereas gains outnumber losses after TNFalpha (802 versus 112).

This comparison evaluates treatment direction across the supplied union peak matrix. The paper's LTR10 enrichment analysis provides the repeat-family-specific context, while this workflow provides DESeq2 region-level testing and nearest-gene annotation.

## Methodological note

DESeq2 median-ratio normalization estimates relative H3K27ac signal and assumes that most tested regions do not share a common directional shift. The supplied matrices have no reported spike-in or IgG/background scaling, so treatment effects are interpreted as relative rather than global changes. Each condition has two biological replicates; the untreated pair shows visible separation in sample-level QC and is retained as observed.

## Citation

Ivancevic A et al. [Endogenous retroviruses mediate transcriptional rewiring in response to oncogenic signaling in colorectal cancer](https://doi.org/10.1126/sciadv.ado1218). *Science Advances*. 2024. Public data and author source files: [GSE186608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186608), [Zenodo 10996183](https://zenodo.org/records/10996183), and pinned source tag `v1.0` / commit `faeef682bebdfb982c70b47a5d7df4ac77f61264`.
