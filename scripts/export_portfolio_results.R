#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))
config <- load_analysis_config(file.path("config", "gse186608_config.R"))
source_dir <- file.path(config$output_dir, "tables")
destination <- file.path("results", "gse186608")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(source_dir, "contrast_summary.csv")
if (!file.exists(summary_path)) stop("Run the GSE186608 analysis before exporting portfolio results.")
summary <- read.csv(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
write.csv(summary, file.path(destination, "contrast_summary.csv"), row.names = FALSE)

snapshots <- lapply(summary$contrast, function(label) {
  path <- file.path(source_dir, paste0(safe_file_stem(label), "_deseq2.csv"))
  result <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  selected <- result[result$test_status == "significant", , drop = FALSE]
  selected <- selected[order(selected$padj, -abs(selected$log2FoldChange), selected$peak_id), , drop = FALSE]
  selected <- head(selected, 20L)
  selected$contrast <- label
  preferred <- c("contrast", "peak_id", "seqnames", "start", "end", "log2FoldChange",
                 "pvalue", "padj", "direction", "annotation", "geneId", "SYMBOL", "distanceToTSS")
  selected[, intersect(preferred, names(selected)), drop = FALSE]
})
write.csv(do.call(rbind, snapshots), file.path(destination, "top_significant_regions.csv"), row.names = FALSE)
cat("Regenerated GSE186608 portfolio summaries.\n")
