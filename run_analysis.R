#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)])
project_root <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
setwd(project_root)
source(file.path("R", "helpers.R"))

if (length(args) && identical(args[[1L]], "--check-dependencies")) {
  config_path <- if (length(args) >= 2L) args[[2L]] else file.path("config", "example_config.R")
  config <- load_analysis_config(config_path)
  dependencies <- check_dependencies(config, include_optional = isTRUE(config$annotation$enabled))
  print(dependencies, row.names = FALSE)
  quit(status = if (any(!dependencies$available)) 1L else 0L)
}

config_path <- if (length(args)) args[[1L]] else file.path("config", "example_config.R")
config <- load_analysis_config(config_path)
config_path_abs <- normalizePath(config_path)
dependencies <- check_dependencies(config, include_optional = isTRUE(config$annotation$enabled))
required_rows <- dependencies$required_for %in% c("core analysis", "core report") |
  (isTRUE(config$annotation$enabled) & dependencies$required_for == "optional genomic annotation")
missing <- dependencies$dependency[required_rows & !dependencies$available]
if (length(missing)) stop("Missing required dependencies: ", paste(missing, collapse = ", "),
                          ". No packages were installed.")

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir_abs <- normalizePath(config$output_dir)
rmarkdown::render(
  input = file.path("report", "cutandrun_downstream.Rmd"),
  output_file = "cutandrun_downstream.html",
  output_dir = output_dir_abs,
  intermediates_dir = output_dir_abs,
  params = list(config_path = config_path_abs),
  knit_root_dir = project_root,
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
