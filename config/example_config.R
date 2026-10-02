analysis_config <- list(
  project_title = "Synthetic CUT&RUN downstream mechanics check",
  data_label = "SYNTHETIC DEMONSTRATION - NO BIOLOGICAL INTERPRETATION",
  counts_path = "example/synthetic/peak_counts.csv",
  metadata_path = "example/synthetic/metadata.csv",
  contrasts_path = "example/synthetic/contrasts.csv",
  output_dir = "output/synthetic",
  coordinate_system = "one_based_closed",
  genome_build = "synthetic_v1",
  design = "~ batch + condition",
  min_total_count = 20L,
  fdr = 0.05,
  pca_top_peaks = 300L,
  heatmap_top_peaks = 30L,
  normalization = list(
    method = "deseq2_median_ratio",
    size_factors_path = NULL
  ),
  annotation = list(
    enabled = FALSE,
    genome_build = "hg38",
    txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    txdb_object = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    orgdb_package = "org.Hs.eg.db",
    tss_region = c(-3000L, 3000L)
  )
)
