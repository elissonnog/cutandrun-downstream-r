required_core_packages <- function() {
  c("rmarkdown", "knitr", "DESeq2", "SummarizedExperiment", "ggplot2")
}

required_annotation_packages <- function(config = NULL) {
  packages <- c("GenomicRanges", "IRanges", "GenomeInfoDb", "AnnotationDbi", "ChIPseeker")
  if (!is.null(config) && nzchar(config$annotation$txdb_package)) {
    packages <- c(packages, config$annotation$txdb_package)
  }
  if (!is.null(config) && nzchar(config$annotation$orgdb_package)) {
    packages <- c(packages, config$annotation$orgdb_package)
  }
  unique(packages)
}

is_pandoc_available <- function() {
  nzchar(Sys.which("pandoc")) ||
    (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available())
}

check_dependencies <- function(config = NULL, include_optional = TRUE) {
  packages <- required_core_packages()
  required_for <- rep("core analysis", length(packages))
  if (isTRUE(include_optional)) {
    optional <- required_annotation_packages(config)
    packages <- c(packages, optional)
    required_for <- c(required_for, rep("optional genomic annotation", length(optional)))
  }
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  data.frame(
    dependency = c(packages, "Pandoc"),
    kind = c(rep("R package", length(packages)), "renderer"),
    required_for = c(required_for, "core report"),
    available = c(unname(available), is_pandoc_available()),
    row.names = NULL
  )
}

load_analysis_config <- function(path) {
  if (!file.exists(path)) stop("Configuration file not found: ", path)
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  if (!exists("analysis_config", envir = env, inherits = FALSE)) {
    stop("Configuration must define an `analysis_config` list.")
  }
  config <- get("analysis_config", envir = env, inherits = FALSE)
  if (!is.list(config)) stop("`analysis_config` must be a list.")
  required <- c(
    "project_title", "data_label", "counts_path", "metadata_path", "contrasts_path",
    "output_dir", "coordinate_system", "genome_build", "design", "min_total_count", "fdr", "pca_top_peaks",
    "heatmap_top_peaks", "normalization", "annotation"
  )
  missing <- setdiff(required, names(config))
  if (length(missing)) stop("Missing configuration fields: ", paste(missing, collapse = ", "))
  if (!is.numeric(config$min_total_count) || length(config$min_total_count) != 1L ||
      is.na(config$min_total_count) || config$min_total_count < 0) {
    stop("`min_total_count` must be one non-negative number.")
  }
  if (!is.numeric(config$fdr) || length(config$fdr) != 1L ||
      is.na(config$fdr) || config$fdr <= 0 || config$fdr >= 1) {
    stop("`fdr` must be one number strictly between 0 and 1.")
  }
  for (field in c("pca_top_peaks", "heatmap_top_peaks")) {
    value <- config[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        value < 2 || value != as.integer(value)) {
      stop("`", field, "` must be one integer greater than or equal to 2.")
    }
  }
  if (!config$coordinate_system %in% c("one_based_closed", "bed_zero_based_half_open")) {
    stop("`coordinate_system` must be `one_based_closed` or `bed_zero_based_half_open`.")
  }
  if (!is.character(config$genome_build) || length(config$genome_build) != 1L ||
      is.na(config$genome_build) || !nzchar(trimws(config$genome_build))) {
    stop("`genome_build` must be one explicit non-blank build identifier.")
  }
  if (!is.list(config$normalization) || !config$normalization$method %in%
      c("deseq2_median_ratio", "provided_deseq2_size_factors")) {
    stop("`normalization$method` must be `deseq2_median_ratio` or `provided_deseq2_size_factors`.")
  }
  if (identical(config$normalization$method, "provided_deseq2_size_factors") &&
      (is.null(config$normalization$size_factors_path) ||
       !nzchar(config$normalization$size_factors_path))) {
    stop("`normalization$size_factors_path` is required for provided size factors.")
  }
  annotation_required <- c(
    "enabled", "genome_build", "txdb_package", "txdb_object",
    "orgdb_package", "tss_region"
  )
  annotation_missing <- setdiff(annotation_required, names(config$annotation))
  if (length(annotation_missing)) {
    stop("Missing annotation configuration fields: ", paste(annotation_missing, collapse = ", "))
  }
  if (!is.numeric(config$annotation$tss_region) || length(config$annotation$tss_region) != 2L ||
      anyNA(config$annotation$tss_region) || config$annotation$tss_region[1L] >= 0 ||
      config$annotation$tss_region[2L] <= 0) {
    stop("`annotation$tss_region` must be two numbers spanning zero, such as c(-3000, 3000).")
  }
  if (isTRUE(config$annotation$enabled) &&
      !identical(config$annotation$genome_build, config$genome_build)) {
    stop("Annotation genome build (`", config$annotation$genome_build,
         "`) must exactly match peak genome build (`", config$genome_build, "`).")
  }
  config
}

read_and_validate_csv_header <- function(path, input_name) {
  header <- read.csv(
    path, header = FALSE, nrows = 1L, check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character"
  )
  header <- unname(unlist(header[1L, ], use.names = FALSE))
  if (anyNA(header) || any(trimws(header) == "")) {
    stop(input_name, " column headers must be non-blank.")
  }
  normalized <- tolower(trimws(header))
  if (anyDuplicated(normalized)) {
    duplicates <- unique(header[duplicated(normalized) | duplicated(normalized, fromLast = TRUE)])
    stop(input_name, " column headers must be unique ignoring case and surrounding spaces; duplicated: ",
         paste(duplicates, collapse = ", "))
  }
  invisible(header)
}

read_peak_counts <- function(path, coordinate_system) {
  if (!file.exists(path)) stop("Peak-count matrix not found: ", path)
  if (!coordinate_system %in% c("one_based_closed", "bed_zero_based_half_open")) {
    stop("Unsupported coordinate system: ", coordinate_system)
  }
  read_and_validate_csv_header(path, "Peak-count matrix")
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  coordinate_columns <- c("peak_id", "seqnames", "start", "end")
  missing <- setdiff(coordinate_columns, names(tab))
  if (length(missing)) {
    stop("Peak-count matrix missing required columns: ", paste(missing, collapse = ", "))
  }
  sample_columns <- setdiff(names(tab), coordinate_columns)
  if (length(sample_columns) < 2L) stop("Peak-count matrix must contain at least two sample columns.")
  if (anyNA(tab$peak_id) || any(trimws(tab$peak_id) == "") || anyDuplicated(tab$peak_id)) {
    stop("`peak_id` values must be non-blank and unique.")
  }
  if (anyNA(tab$seqnames) || any(trimws(tab$seqnames) == "")) {
    stop("`seqnames` values must be non-blank.")
  }
  if (!is.numeric(tab$start) || !is.numeric(tab$end) || anyNA(tab[c("start", "end")]) ||
      any(!is.finite(as.matrix(tab[c("start", "end")]))) ||
      any(tab$start != round(tab$start)) || any(tab$end != round(tab$end))) {
    stop("Peak coordinates must be finite integers.")
  }
  if (identical(coordinate_system, "one_based_closed") &&
      (any(tab$start < 1) || any(tab$end < tab$start))) {
    stop("One-based closed coordinates require `start >= 1` and `end >= start`.")
  }
  if (identical(coordinate_system, "bed_zero_based_half_open") &&
      (any(tab$start < 0) || any(tab$end <= tab$start))) {
    stop("BED zero-based half-open coordinates require `start >= 0` and `end > start`.")
  }
  non_numeric <- sample_columns[!vapply(tab[sample_columns], is.numeric, logical(1))]
  if (length(non_numeric)) stop("Non-numeric count columns: ", paste(non_numeric, collapse = ", "))
  counts <- as.matrix(tab[sample_columns])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- tab$peak_id
  if (anyNA(counts) || any(!is.finite(counts))) stop("Counts must be finite and non-missing.")
  if (any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
    stop("Counts must be non-negative integers.")
  }
  peaks <- tab[coordinate_columns]
  if (identical(coordinate_system, "bed_zero_based_half_open")) {
    peaks$start <- peaks$start + 1L
  }
  coordinate_keys <- paste(peaks$seqnames, peaks$start, peaks$end, sep = ":")
  if (anyDuplicated(coordinate_keys)) {
    stop("Peak coordinate intervals must be unique after coordinate conversion.")
  }
  rownames(peaks) <- peaks$peak_id
  attr(peaks, "input_coordinate_system") <- coordinate_system
  attr(peaks, "analysis_coordinate_system") <- "one_based_closed"
  list(counts = round(counts), peaks = peaks)
}

read_sample_metadata <- function(path) {
  if (!file.exists(path)) stop("Metadata not found: ", path)
  read_and_validate_csv_header(path, "Metadata")
  metadata <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"sample_id" %in% names(metadata)) stop("Metadata must contain a `sample_id` column.")
  if (anyNA(metadata$sample_id) || any(trimws(metadata$sample_id) == "") ||
      anyDuplicated(metadata$sample_id)) {
    stop("`sample_id` values must be non-blank and unique.")
  }
  rownames(metadata) <- metadata$sample_id
  metadata$sample_id <- NULL
  metadata
}

safe_file_stem <- function(x) {
  stem <- gsub("[^A-Za-z0-9._-]+", "_", x)
  stem <- gsub("^_+|_+$", "", stem)
  ifelse(nchar(stem), stem, "contrast")
}

read_contrasts <- function(path) {
  if (!file.exists(path)) stop("Contrast table not found: ", path)
  read_and_validate_csv_header(path, "Contrast table")
  contrasts <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("label", "factor", "numerator", "denominator")
  missing <- setdiff(required, names(contrasts))
  if (length(missing)) stop("Contrast table missing columns: ", paste(missing, collapse = ", "))
  contrasts <- contrasts[required]
  if (!nrow(contrasts)) stop("Contrast table must contain at least one comparison.")
  if (anyNA(contrasts) || any(trimws(as.matrix(contrasts)) == "") || anyDuplicated(contrasts$label)) {
    stop("Contrast values must be complete and labels must be unique.")
  }
  if (any(contrasts$numerator == contrasts$denominator)) {
    labels <- contrasts$label[contrasts$numerator == contrasts$denominator]
    stop("Contrast numerator and denominator must differ; invalid labels: ", paste(labels, collapse = ", "))
  }
  stems <- tolower(safe_file_stem(contrasts$label))
  if (anyDuplicated(stems)) {
    collision_rows <- duplicated(stems) | duplicated(stems, fromLast = TRUE)
    stop("Contrast labels produce case-insensitive colliding output filenames after sanitization: ",
         paste(contrasts$label[collision_rows], collapse = ", "))
  }
  contrasts
}

validate_analysis_inputs <- function(counts, metadata, contrasts, design) {
  if (!identical(colnames(counts), rownames(metadata))) {
    missing_metadata <- setdiff(colnames(counts), rownames(metadata))
    extra_metadata <- setdiff(rownames(metadata), colnames(counts))
    if (length(missing_metadata) || length(extra_metadata)) {
      stop("Count and metadata sample IDs differ. Missing metadata: ",
           paste(missing_metadata, collapse = ", "), "; extra metadata: ",
           paste(extra_metadata, collapse = ", "))
    }
    metadata <- metadata[colnames(counts), , drop = FALSE]
  }
  design_formula <- tryCatch(as.formula(design), error = function(e) NULL)
  if (is.null(design_formula)) stop("Invalid design formula: ", design)
  design_terms <- stats::terms(design_formula)
  if (attr(design_terms, "response") != 0L) stop("The design must not include a response variable.")
  term_labels <- attr(design_terms, "term.labels")
  if (!length(term_labels) || any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", term_labels))) {
    stop("Only additive categorical design terms are supported; interactions, transformations, and continuous covariates are not supported.")
  }
  design_variables <- all.vars(design_formula)
  absent <- setdiff(design_variables, names(metadata))
  if (length(absent)) stop("Design variables absent from metadata: ", paste(absent, collapse = ", "))
  if (anyNA(metadata[design_variables])) stop("Design variables must not contain missing values.")
  non_categorical <- design_variables[!vapply(
    metadata[design_variables],
    function(x) is.character(x) || is.factor(x) || is.logical(x),
    logical(1)
  )]
  if (length(non_categorical)) {
    stop("Design variables must be categorical for this workflow: ", paste(non_categorical, collapse = ", "))
  }
  for (column in design_variables) {
    metadata[[column]] <- factor(metadata[[column]])
    if (nlevels(metadata[[column]]) < 2L) stop("Design variable must contain at least two levels: ", column)
  }
  for (column in unique(contrasts$factor)) {
    if (!column %in% names(metadata)) stop("Contrast factor absent from metadata: ", column)
    if (!column %in% design_variables) stop("Contrast factor absent from design formula: ", column)
    levels_present <- levels(metadata[[column]])
    rows <- contrasts$factor == column
    requested <- unique(c(contrasts$numerator[rows], contrasts$denominator[rows]))
    absent_levels <- setdiff(requested, levels_present)
    if (length(absent_levels)) {
      stop("Contrast levels absent from `", column, "`: ", paste(absent_levels, collapse = ", "))
    }
  }
  model <- stats::model.matrix(design_formula, metadata)
  if (qr(model)$rank < ncol(model)) stop("Design matrix is not full rank; revise confounded categorical terms.")
  metadata
}

read_size_factors <- function(path, sample_ids) {
  if (!file.exists(path)) stop("Size-factor file not found: ", path)
  read_and_validate_csv_header(path, "Size-factor table")
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(tab), c("sample_id", "deseq2_size_factor"))) {
    stop("Size-factor table must contain exactly `sample_id,deseq2_size_factor` columns.")
  }
  if (anyNA(tab) || anyDuplicated(tab$sample_id) || !is.numeric(tab$deseq2_size_factor) ||
      any(!is.finite(tab$deseq2_size_factor)) || any(tab$deseq2_size_factor <= 0)) {
    stop("Provided size factors must be finite, positive, non-missing, with unique sample IDs.")
  }
  if (!setequal(tab$sample_id, sample_ids)) stop("Size-factor and count-matrix sample IDs differ.")
  factors <- tab$deseq2_size_factor[match(sample_ids, tab$sample_id)]
  names(factors) <- sample_ids
  factors / exp(mean(log(factors)))
}

prepare_output_directories <- function(output_dir) {
  paths <- c(root = output_dir, tables = file.path(output_dir, "tables"),
             figures = file.path(output_dir, "figures"))
  for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  paths
}

prepare_deseq_dataset <- function(counts, metadata, design, min_total_count,
                                  normalization, project_root = getwd()) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Package `DESeq2` is required.")
  keep <- rowSums(counts) >= as.numeric(min_total_count)
  if (sum(keep) < 2L) stop("Fewer than two peaks remain after low-count filtering.")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts[keep, , drop = FALSE], colData = metadata, design = as.formula(design)
  )
  if (identical(normalization$method, "provided_deseq2_size_factors")) {
    path <- normalization$size_factors_path
    if (!grepl("^/", path)) path <- file.path(project_root, path)
    DESeq2::sizeFactors(dds) <- read_size_factors(path, colnames(dds))
  }
  dds
}

fit_deseq_model <- function(dds) DESeq2::DESeq(dds, quiet = TRUE)

transform_for_visualization <- function(dds) {
  DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
}

run_deseq_contrast <- function(dds, contrast_row, fdr) {
  result <- DESeq2::results(
    dds,
    contrast = c(contrast_row$factor, contrast_row$numerator, contrast_row$denominator),
    alpha = fdr
  )
  tab <- as.data.frame(result)
  tab$peak_id <- rownames(tab)
  tab$significant <- !is.na(tab$padj) & tab$padj < fdr
  tab$direction <- ifelse(tab$significant & tab$log2FoldChange > 0, "increased",
                          ifelse(tab$significant & tab$log2FoldChange < 0,
                                 "decreased", "not_significant"))
  tab[, c("peak_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue",
          "padj", "significant", "direction")]
}

make_pca_plot <- function(transformed, metadata, top_peaks = 500L, plot_title = "Peak-count PCA") {
  matrix <- SummarizedExperiment::assay(transformed)
  peak_variance <- apply(matrix, 1L, stats::var)
  selected <- names(sort(peak_variance, decreasing = TRUE))[seq_len(min(length(peak_variance), top_peaks))]
  pca <- stats::prcomp(t(matrix[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  data <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1L], PC2 = pca$x[, 2L],
                     metadata[rownames(pca$x), , drop = FALSE], check.names = FALSE)
  color_variable <- if ("condition" %in% names(data)) "condition" else names(metadata)[1L]
  data$display_label <- if ("plot_label" %in% names(data)) data$plot_label else data$sample_id
  ggplot2::ggplot(data, ggplot2::aes(x = PC1, y = PC2, color = .data[[color_variable]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = display_label), vjust = -0.8, size = 3.4,
                       check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::labs(title = plot_title, x = sprintf("PC1 (%.1f%%)", variance[1L]),
                  y = sprintf("PC2 (%.1f%%)", variance[2L]), color = color_variable) +
    ggplot2::theme_minimal(base_size = 11)
}

make_library_qc_plot <- function(counts, metadata, plot_title = "Peak-count library QC") {
  data <- data.frame(
    sample_id = colnames(counts), total_counts = colSums(counts),
    detected_peaks = colSums(counts > 0), metadata[colnames(counts), , drop = FALSE],
    check.names = FALSE
  )
  color_variable <- if ("condition" %in% names(data)) "condition" else names(metadata)[1L]
  data$display_label <- if ("plot_label" %in% names(data)) data$plot_label else data$sample_id
  label_offsets <- c(-1.0, 1.4, -1.0, 1.4, 2.2, -2.0)
  data$label_vjust <- rep(label_offsets, length.out = nrow(data))
  ggplot2::ggplot(data, ggplot2::aes(x = total_counts, y = detected_peaks,
                                    color = .data[[color_variable]], label = display_label)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(vjust = label_vjust), size = 3.4,
                       check_overlap = FALSE, show.legend = FALSE) +
    ggplot2::scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::labs(title = plot_title, x = "Total counts in consensus peaks",
                  y = "Peaks with non-zero counts", color = color_variable) +
    ggplot2::theme_minimal(base_size = 11)
}

make_volcano_plot <- function(result_table, label, fdr, data_label = NULL) {
  plot_data <- result_table
  plot_data$minus_log10_padj <- -log10(pmax(plot_data$padj, .Machine$double.xmin))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = log2FoldChange, y = minus_log10_padj,
                                         color = direction)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(increased = "#B2182B", decreased = "#2166AC",
                                           not_significant = "#BDBDBD")) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed", linewidth = 0.4) +
    ggplot2::labs(title = if (is.null(data_label)) label else paste(data_label, label, sep = ": "),
                  x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}

make_ma_plot <- function(result_table, label, fdr, data_label = NULL) {
  plot_data <- result_table[result_table$baseMean > 0 & is.finite(result_table$log2FoldChange), ]
  ggplot2::ggplot(plot_data, ggplot2::aes(x = baseMean, y = log2FoldChange, color = direction)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(increased = "#B2182B", decreased = "#2166AC",
                                           not_significant = "#BDBDBD")) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    ggplot2::labs(
      title = if (is.null(data_label)) label else paste(data_label, label, sep = ": "),
      subtitle = paste0("Unshrunken DESeq2 estimates; significance is padj < ", fdr),
      x = "Mean normalized count", y = "log2 fold change", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

save_ggplot <- function(plot, path, width = 7, height = 5) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = 150)
  invisible(path)
}

save_correlation_heatmap <- function(transformed, path, plot_title = "Sample correlation", metadata = NULL) {
  correlation <- stats::cor(SummarizedExperiment::assay(transformed), method = "pearson")
  if (!is.null(metadata) && "plot_label" %in% names(metadata)) {
    labels <- as.character(metadata[colnames(correlation), "plot_label"])
    colnames(correlation) <- labels
    rownames(correlation) <- labels
  }
  grDevices::png(path, width = 1400, height = 1200, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  stats::heatmap(correlation, symm = TRUE, margins = c(9, 9), main = plot_title)
  invisible(path)
}

rank_significant_peaks <- function(result_table, fdr, max_peaks) {
  significant <- result_table[!is.na(result_table$padj) & result_table$padj < fdr, , drop = FALSE]
  significant <- significant[order(significant$padj, -abs(significant$log2FoldChange),
                                   significant$peak_id), , drop = FALSE]
  head(significant, max_peaks)
}

save_significant_peak_heatmap <- function(transformed, metadata, path, result_table, fdr,
                                          max_peaks = 30L,
                                          plot_title = "Significant peak replicate heatmap") {
  matrix <- SummarizedExperiment::assay(transformed)
  all_significant <- rank_significant_peaks(result_table, fdr, nrow(result_table))
  selected <- head(intersect(all_significant$peak_id, rownames(matrix)), max_peaks)
  total_significant <- nrow(all_significant)
  grDevices::png(path, width = 1600, height = 1400, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (length(selected) < 2L) {
    graphics::plot.new()
    graphics::title(main = plot_title)
    graphics::text(0.5, 0.5, paste("Fewer than two peaks meet padj <", fdr))
    return(invisible(list(path = path, selected_n = length(selected), significant_n = total_significant)))
  }
  centered <- sweep(matrix[selected, , drop = FALSE], 1L,
                    rowMeans(matrix[selected, , drop = FALSE]), FUN = "-")
  row_order <- stats::hclust(stats::dist(centered))$order
  centered <- centered[row_order, , drop = FALSE]
  condition_variable <- if ("condition" %in% names(metadata)) "condition" else names(metadata)[1L]
  conditions <- as.character(metadata[colnames(centered), condition_variable])
  condition_levels <- unique(conditions)
  condition_palette <- grDevices::hcl.colors(length(condition_levels), "Dark 3")
  names(condition_palette) <- condition_levels
  sample_colors <- unname(condition_palette[conditions])
  limit <- max(abs(centered))
  heat_colors <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101L)

  graphics::layout(matrix(c(1, 2, 3), ncol = 1L), heights = c(0.35, 5, 0.75))
  graphics::par(oma = c(0, 0, 3.3, 0), mar = c(0, 5, 1, 13))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, ncol(centered) + 0.5), ylim = c(0, 1))
  for (index in seq_len(ncol(centered))) {
    graphics::rect(index - 0.5, 0, index + 0.5, 1,
                   col = sample_colors[index], border = NA)
  }
  graphics::mtext("Condition", side = 2, line = 2.5, cex = 0.8)

  graphics::par(mar = c(6, 5, 0.5, 13))
  graphics::image(
    x = seq_len(ncol(centered)), y = seq_len(nrow(centered)), z = t(centered),
    col = heat_colors, zlim = c(-limit, limit), axes = FALSE, xlab = "", ylab = "",
    useRaster = TRUE
  )
  sample_labels <- if ("plot_label" %in% names(metadata)) {
    as.character(metadata[colnames(centered), "plot_label"])
  } else {
    colnames(centered)
  }
  graphics::axis(1, at = seq_len(ncol(centered)), labels = sample_labels, las = 2, cex.axis = 0.8)
  graphics::axis(4, at = seq_len(nrow(centered)), labels = rownames(centered),
                 las = 1, cex.axis = 0.55, tick = FALSE)

  graphics::par(mar = c(1, 5, 0, 13))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  key_x <- seq(0.05, 0.43, length.out = length(heat_colors) + 1L)
  for (index in seq_along(heat_colors)) {
    graphics::rect(key_x[index], 0.38, key_x[index + 1L], 0.64,
                   col = heat_colors[index], border = NA)
  }
  graphics::text(c(0.05, 0.24, 0.43), 0.2,
                 labels = format(c(-limit, 0, limit), digits = 2), cex = 0.75)
  graphics::text(0.24, 0.86, "Row-centered VST", cex = 0.8)
  graphics::legend("right", legend = condition_levels, fill = condition_palette,
                   title = condition_variable, bty = "n", horiz = TRUE, cex = 0.75)
  graphics::mtext(plot_title, outer = TRUE, side = 3, line = 1.6, cex = 1.15, font = 2)
  graphics::mtext(
    sprintf("Top %d of %d peaks with padj < %s; rows clustered, samples shown in input order",
            length(selected), total_significant, format(fdr)),
    outer = TRUE, side = 3, line = 0.25, cex = 0.78
  )
  invisible(list(path = path, selected_n = length(selected), significant_n = total_significant))
}

annotate_peaks_optional <- function(peaks, config) {
  packages <- required_annotation_packages(config)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Optional annotation requested but packages are unavailable: ", paste(missing, collapse = ", "))
  }
  txdb <- getExportedValue(config$annotation$txdb_package, config$annotation$txdb_object)
  txdb_builds <- unique(stats::na.omit(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(txdb))))
  if (!length(txdb_builds)) {
    stop("Selected TxDb does not declare a genome build; compatibility cannot be verified.")
  }
  if (!identical(config$annotation$genome_build, config$genome_build) ||
      !config$genome_build %in% txdb_builds) {
    stop("Selected TxDb build(s) [", paste(txdb_builds, collapse = ", "),
         "] do not match declared peak genome build `", config$genome_build, "`.")
  }
  incompatible <- setdiff(unique(peaks$seqnames), GenomeInfoDb::seqlevels(txdb))
  if (length(incompatible)) {
    stop("Peak seqnames absent from the selected TxDb: ",
         paste(head(incompatible, 10L), collapse = ", "),
         if (length(incompatible) > 10L) " ..." else "")
  }
  ranges <- GenomicRanges::GRanges(
    seqnames = peaks$seqnames,
    ranges = IRanges::IRanges(start = peaks$start, end = peaks$end),
    peak_id = peaks$peak_id
  )
  names(ranges) <- peaks$peak_id
  annotated <- ChIPseeker::annotatePeak(
    ranges, TxDb = txdb, annoDb = config$annotation$orgdb_package,
    tssRegion = as.integer(config$annotation$tss_region), verbose = FALSE
  )
  annotation <- as.data.frame(annotated)
  annotation$peak_id <- NULL
  joined <- merge(peaks, annotation, by = c("seqnames", "start", "end"),
                  all.x = TRUE, sort = FALSE)
  joined <- joined[match(peaks$peak_id, joined$peak_id), , drop = FALSE]
  if (!identical(joined$peak_id, peaks$peak_id)) {
    stop("Annotation join did not preserve unique peak identifiers and order.")
  }
  joined
}
