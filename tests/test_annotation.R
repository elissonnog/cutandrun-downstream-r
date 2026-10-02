#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

assert_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NA_character_ },
                      error = function(error) conditionMessage(error))
  if (is.na(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop("Expected error containing `", pattern, "`; received: ", message, call. = FALSE)
  }
}

config <- load_analysis_config(file.path("config", "annotation_test_config.R"))
dependencies <- check_dependencies(config, include_optional = TRUE)
missing <- dependencies$dependency[dependencies$required_for == "optional genomic annotation" &
                                     !dependencies$available]
if (length(missing)) stop("Annotation test requires: ", paste(missing, collapse = ", "))

fixture <- read_peak_counts(config$counts_path, config$coordinate_system)
annotation <- annotate_peaks_optional(fixture$peaks, config)

assert_true(nrow(annotation) == nrow(fixture$peaks), "Annotation changed the peak count.")
assert_true(identical(annotation$peak_id, fixture$peaks$peak_id),
            "Annotation changed peak IDs or order.")
assert_true(all(c("annotation", "geneId", "SYMBOL") %in% names(annotation)),
            "Expected ChIPseeker/OrgDb annotation columns are missing.")
assert_true(sum(!is.na(annotation$SYMBOL)) > 0L, "No reference gene symbols were assigned.")

bad_seqnames <- fixture$peaks[1L, , drop = FALSE]
bad_seqnames$seqnames <- "chrSynthetic"
assert_error(annotate_peaks_optional(bad_seqnames, config),
             "Peak seqnames absent from the selected TxDb")

bed_path <- tempfile(fileext = ".csv")
bed_counts <- data.frame(
  peak_id = fixture$peaks$peak_id[1:3],
  seqnames = fixture$peaks$seqnames[1:3],
  start = fixture$peaks$start[1:3] - 1L,
  end = fixture$peaks$end[1:3],
  SampleA = c(10L, 20L, 30L),
  SampleB = c(11L, 21L, 31L)
)
write.csv(bed_counts, bed_path, row.names = FALSE, quote = FALSE)
bed_fixture <- read_peak_counts(bed_path, "bed_zero_based_half_open")
assert_true(identical(bed_fixture$peaks$start, fixture$peaks$start[1:3]),
            "BED conversion did not recover the reference one-based starts.")
bed_annotation <- annotate_peaks_optional(bed_fixture$peaks, config)
assert_true(identical(bed_annotation$peak_id, bed_fixture$peaks$peak_id),
            "BED-derived annotation did not preserve peak IDs.")

txdb <- getExportedValue(config$annotation$txdb_package, config$annotation$txdb_object)
txdb_builds <- unique(stats::na.omit(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(txdb))))
assert_true(identical(txdb_builds, "hg38"), "TxDb does not report the expected hg38 build.")

cat("Annotation mechanics test passed for", nrow(annotation), "reference-grounded peaks.\n")
