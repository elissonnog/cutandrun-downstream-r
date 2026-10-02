#!/usr/bin/env Rscript

# Download and reshape the pinned public GSE186608 H3K27ac CUT&RUN example.
# Source archive: https://zenodo.org/records/10996183 (CC BY 4.0)
# Modifications here: checksum verification, exact BED/count mapping validation,
# schema conversion, and explicit sample/contrast metadata for this workflow.

if (!requireNamespace("openssl", quietly = TRUE)) {
  stop("Package `openssl` is required for SHA-256 verification.")
}

manifest <- read.csv(file.path("config", "gse186608_sources.csv"), stringsAsFactors = FALSE)
source_dir <- file.path("data", "gse186608", "source")
input_dir <- file.path("data", "gse186608", "analysis_input")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  unclass(unname(as.character(openssl::sha256(connection))))
}

for (index in seq_len(nrow(manifest))) {
  destination <- file.path(source_dir, manifest$file[index])
  if (!file.exists(destination)) {
    utils::download.file(manifest$url[index], destination, mode = "wb", quiet = FALSE)
  }
  observed_bytes <- file.info(destination)$size
  observed_sha256 <- sha256_file(destination)
  if (!identical(as.numeric(observed_bytes), as.numeric(manifest$bytes[index]))) {
    stop("Byte-size mismatch for ", manifest$file[index], ": expected ",
         manifest$bytes[index], ", observed ", observed_bytes)
  }
  if (!identical(observed_sha256, unname(manifest$sha256[index]))) {
    stop("SHA-256 mismatch for ", manifest$file[index], ": expected ",
         manifest$sha256[index], ", observed ", observed_sha256)
  }
}

counts_path <- file.path(source_dir, "raw_counts_H3K27ac_CnR.tab")
bed_path <- file.path(source_dir, "ALL_K27ac_CnR_top20k_peaks.bed")
counts <- read.delim(counts_path, header = TRUE, row.names = 1L,
                     check.names = FALSE, stringsAsFactors = FALSE)
bed <- read.delim(bed_path, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("seqnames", "start", "end", "peak_id"))

expected_samples <- c(
  "UT_K27ac_rep1", "UT_K27ac_rep2", "Cobi_K27ac_rep1",
  "Cobi_K27ac_rep2", "TNFa_K27ac_rep1", "TNFa_K27ac_rep2"
)
coordinate_ids <- paste0(bed$seqnames, ":", bed$start, "-", bed$end)
if (!identical(dim(counts), c(38649L, 6L)) || nrow(bed) != 38649L) {
  stop("Pinned input dimensions changed; expected 38,649 regions by six samples.")
}
if (!identical(colnames(counts), expected_samples)) stop("Pinned count headers changed.")
if (!identical(rownames(counts), coordinate_ids)) {
  stop("Count row coordinate IDs do not exactly match BED rows in order.")
}
if (anyDuplicated(coordinate_ids) || anyDuplicated(bed$peak_id)) {
  stop("Pinned BED contains duplicated intervals or peak IDs.")
}
count_matrix <- as.matrix(counts)
if (anyNA(count_matrix) || any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
    any(count_matrix != round(count_matrix))) {
  stop("Pinned counts must be finite, non-missing, non-negative integers.")
}
if (any(bed$start < 0) || any(bed$end <= bed$start)) {
  stop("Pinned BED coordinates must be zero-based half-open with end > start.")
}

peak_counts <- data.frame(bed, counts, check.names = FALSE)
metadata <- data.frame(
  sample_id = expected_samples,
  condition = c("untreated", "untreated", "cobimetinib", "cobimetinib", "TNFalpha", "TNFalpha"),
  plot_label = c("UT1", "UT2", "Cobi1", "Cobi2", "TNFa1", "TNFa2"),
  geo_accession = c("GSM5657497", "GSM5657498", "GSM5657501", "GSM5657502", "GSM5657505", "GSM5657506"),
  source_note = c(
    "treatment replicate 1", "treatment replicate 2; source reports resequenced/merged",
    "treatment replicate 1", "treatment replicate 2",
    "treatment replicate 1", "treatment replicate 2"
  ),
  stringsAsFactors = FALSE
)
contrasts <- data.frame(
  label = c("cobimetinib_vs_untreated", "TNFalpha_vs_untreated"),
  factor = "condition",
  numerator = c("cobimetinib", "TNFalpha"),
  denominator = "untreated",
  stringsAsFactors = FALSE
)

write.csv(peak_counts, file.path(input_dir, "peak_counts.csv"), row.names = FALSE, quote = FALSE)
write.csv(metadata, file.path(input_dir, "metadata.csv"), row.names = FALSE, quote = TRUE)
write.csv(contrasts, file.path(input_dir, "contrasts.csv"), row.names = FALSE, quote = FALSE)

cat("Prepared GSE186608 H3K27ac inputs:\n")
cat("  regions:", nrow(peak_counts), "(including", sum(rowSums(count_matrix) == 0), "all-zero regions)\n")
cat("  samples:", ncol(count_matrix), "\n")
cat("  exact BED/count mapping: TRUE\n")
cat("  source checksums: verified\n")
