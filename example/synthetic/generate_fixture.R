#!/usr/bin/env Rscript

set.seed(20261002)

output_dir <- file.path("example", "synthetic")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_peaks <- 800L
sample_id <- sprintf("Demo%02d", seq_len(8L))
condition <- rep(c("reference", "perturbed"), each = 4L)
batch <- rep(c("batch_a", "batch_b"), times = 4L)

# An overdispersed negative-binomial mechanics fixture. Coordinates are synthetic
# and intentionally use chrSynthetic, so they cannot be mistaken for a genome build.
base_mean <- exp(rnorm(n_peaks, log(70), 0.9))
dispersion <- runif(n_peaks, 0.12, 0.32)
fold_change <- rep(1, n_peaks)
fold_change[1:80] <- 3.0
fold_change[81:160] <- 1 / 3.0
library_factor <- c(0.82, 1.05, 0.93, 1.18, 0.88, 1.12, 0.98, 1.25)
batch_factor <- ifelse(batch == "batch_b", 1.10, 1.00)

counts <- vapply(seq_along(sample_id), function(index) {
  condition_factor <- if (condition[index] == "perturbed") fold_change else 1
  mu <- base_mean * condition_factor * library_factor[index] * batch_factor[index]
  stats::rnbinom(n_peaks, mu = mu, size = 1 / dispersion)
}, numeric(n_peaks))
colnames(counts) <- sample_id

starts <- 1000L + (seq_len(n_peaks) - 1L) * 500L
count_table <- data.frame(
  peak_id = sprintf("SyntheticPeak%04d", seq_len(n_peaks)),
  seqnames = "chrSynthetic",
  start = starts,
  end = starts + 249L,
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

cat("Wrote deterministic synthetic CUT&RUN mechanics fixture to", output_dir, "\n")
