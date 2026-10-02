#!/usr/bin/env Rscript

config <- new.env(parent = baseenv())
sys.source(file.path("config", "example_config.R"), envir = config)
bad <- config$analysis_config
bad$annotation$enabled <- TRUE
bad$annotation$genome_build <- bad$genome_build
bad$annotation$txdb_package <- "DefinitelyMissingPackageForPreflight"
bad$annotation$txdb_object <- "DefinitelyMissingPackageForPreflight"
bad$annotation$orgdb_package <- "DefinitelyMissingPackageForPreflight"
path <- tempfile(fileext = ".R")
writeLines(c("analysis_config <-", capture.output(dput(bad))), path)
rscript <- file.path(R.home("bin"), "Rscript")
status_ok <- system2(rscript, c("run_analysis.R", "--check-dependencies", "config/example_config.R"),
                     stdout = FALSE, stderr = FALSE)
status_bad <- system2(rscript, c("run_analysis.R", "--check-dependencies", path),
                      stdout = FALSE, stderr = FALSE)
if (!identical(status_ok, 0L)) stop("Dependency preflight rejected the available synthetic configuration.")
if (identical(status_bad, 0L)) stop("Dependency preflight did not fail for a missing configured package.")
cat("Dependency preflight exit-status test passed.\n")
