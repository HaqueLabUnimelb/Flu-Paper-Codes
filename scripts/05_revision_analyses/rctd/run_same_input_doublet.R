#!/usr/bin/env Rscript
# Purpose: rerun a fitted day-4 RCTD.replicates object in doublet mode.
# Input: FLU_RCTD_DAY4_FULL_RDS.
# Output: FLU_RCTD_DAY4_DOUBLET_RDS.

suppressPackageStartupMessages(library(spacexr))

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Unable to determine the script location. Run this file with Rscript.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
repo_dir <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)

env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) {
    stop("Required path is empty: ", name)
  }
  value
}

input_file <- env_path(
  "FLU_RCTD_DAY4_FULL_RDS",
  file.path(repo_dir, "outputs", "rctd", "Day4mouseRCTDoutmin30_multi.rds")
)
output_file <- env_path(
  "FLU_RCTD_DAY4_DOUBLET_RDS",
  file.path(repo_dir, "outputs", "revision", "rctd", "Day4mouseRCTDoutmin30_doublet_same_input.rds")
)

if (!file.exists(input_file)) {
  stop("Full-mode RCTD input was not found: ", input_file)
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

message("R: ", R.version.string)
message("spacexr: ", as.character(packageVersion("spacexr")))
message("Loading same-input full object: ", input_file)

rctd_doublet <- readRDS(input_file)
if (!inherits(rctd_doublet, "RCTD.replicates")) {
  stop("Expected RCTD.replicates, got: ", paste(class(rctd_doublet), collapse = ", "))
}

message("Re-running identical RCTD inputs with doublet_mode='doublet'")
rctd_doublet <- run.RCTD.replicates(rctd_doublet, doublet_mode = "doublet")
saveRDS(rctd_doublet, output_file)
message("Saved same-input doublet result: ", output_file)
