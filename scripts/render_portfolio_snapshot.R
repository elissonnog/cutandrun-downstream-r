#!/usr/bin/env Rscript

config_path <- normalizePath(file.path("config", "gse186608_config.R"))
rmarkdown::render(
  input = file.path("report", "cutandrun_downstream.Rmd"),
  output_file = "gse186608_h3k27ac_downstream.html",
  output_dir = normalizePath("report"),
  params = list(config_path = config_path, portfolio_snapshot = TRUE),
  knit_root_dir = getwd(), envir = new.env(parent = globalenv()), quiet = FALSE
)
