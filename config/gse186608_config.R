analysis_config <- list(
  project_title = "GSE186608 HCT116 H3K27ac CUT&RUN downstream reanalysis",
  data_label = "PUBLIC DATA REANALYSIS - GSE186608 HCT116 H3K27ac",
  figure_label = "GSE186608 HCT116 H3K27ac",
  counts_path = "data/gse186608/analysis_input/peak_counts.csv",
  metadata_path = "data/gse186608/analysis_input/metadata.csv",
  contrasts_path = "data/gse186608/analysis_input/contrasts.csv",
  output_dir = "output/gse186608_h3k27ac",
  coordinate_system = "bed_zero_based_half_open",
  genome_build = "hg38",
  design = "~ condition",
  min_total_count = 0L,
  fdr = 0.05,
  pca_top_peaks = 500L,
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
  ),
  provenance = list(
    accession = "GSE186608",
    citation = "Ivancevic et al., Science Advances (2024)",
    doi = "10.1126/sciadv.ado1218",
    archive = "https://zenodo.org/records/10996183",
    source_tag = "v1.0",
    source_commit = "faeef682bebdfb982c70b47a5d7df4ac77f61264",
    source_license = "CC BY 4.0",
    counting_unit = "bedtools multicov aligned-read overlaps",
    region_set = "union of per-sample top-20k peak sets",
    inference_limit = "n=2 treatment replicates per condition; relative signal only"
  )
)
