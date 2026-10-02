analysis_config <- list(
  project_title = "Synthetic CUT&RUN hg38 annotation mechanics check",
  data_label = "SYNTHETIC COUNTS + PUBLIC hg38 REFERENCE COORDINATES - NO BIOLOGICAL INTERPRETATION",
  counts_path = "example/synthetic_hg38_annotation/peak_counts.csv",
  metadata_path = "example/synthetic_hg38_annotation/metadata.csv",
  contrasts_path = "example/synthetic_hg38_annotation/contrasts.csv",
  output_dir = "output/synthetic_hg38_annotation",
  coordinate_system = "one_based_closed",
  genome_build = "hg38",
  design = "~ batch + condition",
  min_total_count = 20L,
  fdr = 0.05,
  pca_top_peaks = 200L,
  heatmap_top_peaks = 30L,
  normalization = list(
    method = "deseq2_median_ratio",
    size_factors_path = NULL
  ),
  annotation = list(
    enabled = TRUE,
    genome_build = "hg38",
    txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    txdb_object = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    orgdb_package = "org.Hs.eg.db",
    tss_region = c(-3000L, 500L)
  )
)
