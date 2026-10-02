#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

assert_error <- function(expression, pattern) {
  message <- tryCatch({
    force(expression)
    NA_character_
  }, error = function(error) conditionMessage(error))
  if (is.na(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop("Expected error containing `", pattern, "`; received: ", message, call. = FALSE)
  }
}

fixture <- read_peak_counts(file.path("example", "synthetic", "peak_counts.csv"), "one_based_closed")
metadata <- read_sample_metadata(file.path("example", "synthetic", "metadata.csv"))
contrasts <- read_contrasts(file.path("example", "synthetic", "contrasts.csv"))
metadata <- validate_analysis_inputs(fixture$counts, metadata, contrasts, "~ batch + condition")

assert_true(identical(dim(fixture$counts), c(800L, 8L)), "Synthetic count dimensions changed.")
assert_true(identical(colnames(fixture$counts), rownames(metadata)), "Sample alignment failed.")
assert_true(all(fixture$counts >= 0 & fixture$counts == round(fixture$counts)), "Fixture counts are invalid.")
assert_true(identical(safe_file_stem("treated vs reference"), "treated_vs_reference"),
            "Filename sanitization failed.")

ranking_fixture <- data.frame(
  peak_id = c("peak_b", "peak_a", "peak_c", "peak_d"),
  padj = c(0.01, 0.01, 0.20, NA_real_),
  log2FoldChange = c(1, -2, 8, 10)
)
ranked_significant <- rank_significant_peaks(ranking_fixture, fdr = 0.05, max_peaks = 30L)
assert_true(identical(ranked_significant$peak_id, c("peak_a", "peak_b")),
            "Significant-peak ranking is not deterministic or included non-significant peaks.")

reordered <- metadata[rev(rownames(metadata)), , drop = FALSE]
realigned <- validate_analysis_inputs(fixture$counts, reordered, contrasts, "~ batch + condition")
assert_true(identical(colnames(fixture$counts), rownames(realigned)), "Metadata was not safely reordered.")

bad_contrasts <- contrasts
bad_contrasts$numerator <- "absent_level"
assert_error(validate_analysis_inputs(fixture$counts, metadata, bad_contrasts, "~ batch + condition"),
             "Contrast levels absent")
assert_error(validate_analysis_inputs(fixture$counts, metadata, contrasts, "~ missing_variable"),
             "Design variables absent")
assert_error(validate_analysis_inputs(fixture$counts, metadata, contrasts, "~ batch"),
             "Contrast factor absent from design formula")
assert_error(validate_analysis_inputs(fixture$counts, metadata, contrasts, "~ batch * condition"),
             "Only additive categorical design terms")

numeric_metadata <- metadata
numeric_metadata$continuous <- seq_len(nrow(numeric_metadata))
assert_error(validate_analysis_inputs(fixture$counts, numeric_metadata, contrasts,
                                      "~ batch + condition + continuous"),
             "Design variables must be categorical")

write_temp <- function(lines) {
  path <- tempfile(fileext = ".csv")
  writeLines(lines, path)
  path
}

assert_error(read_contrasts(write_temp(c(
  "label,factor,numerator,denominator",
  "A/B,condition,perturbed,reference",
  "a b,condition,perturbed,reference"
))), "case-insensitive colliding output filenames")

assert_error(read_contrasts(write_temp(c(
  "label,factor,numerator,denominator",
  "same,condition,reference,reference"
))), "numerator and denominator must differ")

assert_error(read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,Demo01,demo01",
  "P1,chrSynthetic,1,10,1,2"
)), "one_based_closed"), "column headers must be unique")

assert_error(read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,,Demo02",
  "P1,chrSynthetic,1,10,1,2"
)), "one_based_closed"), "column headers must be non-blank")

assert_error(read_sample_metadata(write_temp(c(
  "sample_id,condition,Condition",
  "Demo01,reference,reference"
))), "column headers must be unique")

assert_error(read_sample_metadata(write_temp(c(
  "sample_id,,batch",
  "Demo01,reference,batch_a"
))), "column headers must be non-blank")

assert_error(read_contrasts(write_temp(c(
  "label,factor,numerator,Numerator",
  "comparison,condition,perturbed,reference"
))), "column headers must be unique")

assert_error(read_contrasts(write_temp(c(
  "label,factor,,denominator",
  "comparison,condition,perturbed,reference"
))), "column headers must be non-blank")

assert_error(read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,Demo01,Demo02",
  "P1,chrSynthetic,1,10,1.5,2"
)), "one_based_closed"), "Counts must be non-negative integers")

assert_error(read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,Demo01,Demo02",
  "P1,chrSynthetic,1,10,-1,2"
)), "one_based_closed"), "Counts must be non-negative integers")

assert_error(read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,Demo01,Demo02",
  "P1,chrSynthetic,10,1,1,2"
)), "one_based_closed"), "One-based closed coordinates")

bed_fixture <- read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,Demo01,Demo02",
  "P1,chrSynthetic,0,10,1,2"
)), "bed_zero_based_half_open")
assert_true(isTRUE(all.equal(bed_fixture$peaks$start, 1)), "BED start was not converted to one-based closed.")
assert_true(isTRUE(all.equal(bed_fixture$peaks$end, 10)), "BED end changed during conversion.")
assert_error(read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,Demo01,Demo02",
  "P1,chrSynthetic,0,0,1,2"
)), "bed_zero_based_half_open"), "BED zero-based half-open coordinates")
assert_error(read_peak_counts(write_temp(c(
  "peak_id,seqnames,start,end,Demo01,Demo02",
  "P1,chrSynthetic,1,10,1,2"
)), "zero_based_closed"), "Unsupported coordinate system")

size_factor_file <- write_temp(c(
  "sample_id,deseq2_size_factor",
  paste(rev(colnames(fixture$counts)), rev(seq(0.8, 1.5, length.out = ncol(fixture$counts))), sep = ",")
))
size_factors <- read_size_factors(size_factor_file, colnames(fixture$counts))
assert_true(all(is.finite(size_factors) & size_factors > 0), "Valid size factors were rejected.")
assert_true(abs(exp(mean(log(size_factors))) - 1) < 1e-10, "Size factors were not centered.")
assert_true(identical(names(size_factors), colnames(fixture$counts)), "Size factors were not sample-aligned.")
expected_centered <- seq(0.8, 1.5, length.out = ncol(fixture$counts))
expected_centered <- expected_centered / exp(mean(log(expected_centered)))
assert_true(isTRUE(all.equal(unname(size_factors), expected_centered)), "Size factors were assigned out of order.")

assert_error(read_size_factors(write_temp(c(
  "sample_id,deseq2_size_factor",
  "Demo01,0",
  "Demo02,1"
)), c("Demo01", "Demo02")), "finite, positive")

config <- load_analysis_config(file.path("config", "example_config.R"))
assert_true(identical(config$design, "~ batch + condition"), "Configuration did not load correctly.")
assert_true(grepl("SYNTHETIC", config$data_label, fixed = TRUE), "Synthetic label is missing.")
assert_true(!isTRUE(config$annotation$enabled), "Synthetic annotation must remain disabled.")
assert_true(identical(config$genome_build, "synthetic_v1"), "Synthetic genome build is not explicit.")

bad_config <- config
bad_config$annotation$enabled <- TRUE
bad_config_path <- tempfile(fileext = ".R")
writeLines(c("analysis_config <-", capture.output(dput(bad_config))), bad_config_path)
assert_error(load_analysis_config(bad_config_path), "must exactly match peak genome build")

cat("All dependency-free validation and fixture tests passed.\n")
