#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

missing <- required_core_packages()[!vapply(required_core_packages(), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Core integration test requires: ", paste(missing, collapse = ", "))

config <- load_analysis_config(file.path("config", "example_config.R"))
fixture <- read_peak_counts(config$counts_path, config$coordinate_system)
metadata <- read_sample_metadata(config$metadata_path)
contrasts <- read_contrasts(config$contrasts_path)
metadata <- validate_analysis_inputs(fixture$counts, metadata, contrasts, config$design)

dds <- prepare_deseq_dataset(
  fixture$counts, metadata, config$design, config$min_total_count,
  config$normalization, project_root = getwd()
)
dds <- fit_deseq_model(dds)
result <- run_deseq_contrast(dds, contrasts[1L, , drop = FALSE], config$fdr,
                             input_counts = fixture$counts)

assert_true(nrow(result) == 800L, "Integration result lost synthetic peaks.")
assert_true(all(is.finite(DESeq2::sizeFactors(dds)) & DESeq2::sizeFactors(dds) > 0),
            "DESeq2 size factors are invalid.")
assert_true(sum(result$significant, na.rm = TRUE) >= 80L,
            "Synthetic planted signal was not recovered; dispersion/model behavior changed.")
assert_true(stats::median(result$log2FoldChange[match(sprintf("SyntheticPeak%04d", 1:80), result$peak_id)],
                          na.rm = TRUE) > 1,
            "Planted increased peaks do not have the expected effect direction.")
assert_true(stats::median(result$log2FoldChange[match(sprintf("SyntheticPeak%04d", 81:160), result$peak_id)],
                          na.rm = TRUE) < -1,
            "Planted decreased peaks do not have the expected effect direction.")
assert_true(all(result$test_status %in% c(
  "all_zero", "low_count_filtered", "unavailable", "independent_filtered",
  "tested_nonsignificant", "significant"
)), "Unexpected result status.")
assert_true(all(is.na(result$significant[result$test_status %in% c(
  "all_zero", "low_count_filtered", "unavailable", "independent_filtered"
)])), "Untested or independently filtered peaks must have significant = NA.")
assert_true(all(!result$significant[result$test_status == "tested_nonsignificant"]),
            "Tested nonsignificant peaks must have significant = FALSE.")

provided <- seq(0.75, 1.45, length.out = ncol(fixture$counts))
factor_path <- tempfile(fileext = ".csv")
write.csv(data.frame(
  sample_id = rev(colnames(fixture$counts)),
  deseq2_size_factor = rev(provided)
), factor_path, row.names = FALSE, quote = FALSE)
provided_normalization <- list(
  method = "provided_deseq2_size_factors",
  size_factors_path = factor_path
)
provided_dds <- prepare_deseq_dataset(
  fixture$counts, metadata, config$design, config$min_total_count,
  provided_normalization, project_root = getwd()
)
expected <- provided / exp(mean(log(provided)))
assert_true(identical(names(DESeq2::sizeFactors(provided_dds)), colnames(fixture$counts)),
            "Assigned DESeq2 factors lost sample ordering.")
assert_true(isTRUE(all.equal(unname(DESeq2::sizeFactors(provided_dds)), expected)),
            "Assigned DESeq2 factors do not match the declared divisors.")

cat("Core DESeq2 integration test passed.\n")
