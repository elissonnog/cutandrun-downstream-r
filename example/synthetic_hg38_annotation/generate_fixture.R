#!/usr/bin/env Rscript

# Public-reference-grounded hg38 coordinates plus synthetic counts. This fixture
# tests annotation mechanics only and is not an experimental or biological dataset.
required <- c("GenomicFeatures", "GenomeInfoDb", "TxDb.Hsapiens.UCSC.hg38.knownGene")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Fixture generation requires: ", paste(missing, collapse = ", "))

set.seed(20261003)
output_dir <- file.path("example", "synthetic_hg38_annotation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

txdb <- getExportedValue(
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
)
genes <- suppressMessages(GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE))
genes <- GenomeInfoDb::keepStandardChromosomes(genes, pruning.mode = "coarse")
gene_table <- as.data.frame(genes)
gene_table <- gene_table[gene_table$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]
gene_table <- gene_table[order(gene_table$seqnames, gene_table$start, gene_table$end), ]

n_peaks <- 240L
gene_table <- gene_table[seq_len(n_peaks), ]
tss <- ifelse(gene_table$strand == "-", gene_table$end, gene_table$start)
peak_start <- pmax(1L, as.integer(tss) - 125L)
peak_end <- as.integer(tss) + 124L

sample_id <- sprintf("Hg%02d", seq_len(8L))
condition <- rep(c("reference", "perturbed"), each = 4L)
batch <- rep(c("batch_a", "batch_b"), times = 4L)
base_mean <- exp(rnorm(n_peaks, log(80), 0.8))
dispersion <- runif(n_peaks, 0.12, 0.28)
fold_change <- rep(1, n_peaks)
fold_change[1:30] <- 2.8
fold_change[31:60] <- 1 / 2.8
library_factor <- c(0.85, 1.08, 0.94, 1.16, 0.90, 1.12, 1.00, 1.22)

counts <- vapply(seq_along(sample_id), function(index) {
  condition_factor <- if (condition[index] == "perturbed") fold_change else 1
  mu <- base_mean * condition_factor * library_factor[index]
  stats::rnbinom(n_peaks, mu = mu, size = 1 / dispersion)
}, numeric(n_peaks))
colnames(counts) <- sample_id

count_table <- data.frame(
  peak_id = sprintf("SyntheticHg38Peak%03d", seq_len(n_peaks)),
  seqnames = as.character(gene_table$seqnames),
  start = peak_start,
  end = peak_end,
  counts,
  check.names = FALSE
)
metadata <- data.frame(sample_id = sample_id, condition = condition, batch = batch)
contrasts <- data.frame(
  label = "perturbed_vs_reference",
  factor = "condition",
  numerator = "perturbed",
  denominator = "reference"
)

write.csv(count_table, file.path(output_dir, "peak_counts.csv"), row.names = FALSE, quote = FALSE)
write.csv(metadata, file.path(output_dir, "metadata.csv"), row.names = FALSE, quote = FALSE)
write.csv(contrasts, file.path(output_dir, "contrasts.csv"), row.names = FALSE, quote = FALSE)
cat("Wrote public-reference-grounded hg38 coordinates with synthetic counts to", output_dir, "\n")
