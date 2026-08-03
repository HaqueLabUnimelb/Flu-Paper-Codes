#!/usr/bin/env Rscript
# Exploratory visualization of cell-type weights from a multi-replicate RCTD result.
# By default the script reads the day-10 result from outputs/rctd. Set
# FLU_RCTD_VIS_INPUT_RDS to select a different result. The original script did
# not save a file; set FLU_RCTD_VIS_OUTPUT to request an optional plot file.

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(tidyverse)
  library(data.table)
  library(mltools)
})

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Unable to determine the script location. Run this file with Rscript.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

input_file <- Sys.getenv(
  "FLU_RCTD_VIS_INPUT_RDS",
  unset = file.path(repo_dir, "outputs", "rctd", "Day10mouseRCTDoutmin30_multi.rds")
)
if (!file.exists(input_file)) {
  stop("RCTD result was not found: ", input_file)
}

normalize <- function(x) {
  x <- (x - base::min(x)) / (max(x) - base::min(x))
  x
}

rctd.res <- readRDS(input_file)
rctd.res <- rctd.res@RCTD.reps
if (length(rctd.res) < 2L) {
  stop("The RCTD result must contain at least two replicates.")
}

res <- list()
for (i in 1:2) {
  map <- as.data.frame(as.matrix(rctd.res[[i]]@results$weights))
  map <- as.data.frame(apply(X = map, MARGIN = 2, FUN = normalize))

  coords <- as.data.frame(rctd.res[[i]]@spatialRNA@coords)
  coords <- coords[rownames(map), ]
  coords$x <- coords$x - mean(coords$x)
  coords$y <- coords$y - mean(coords$y)

  map$x <- coords$x
  map$y <- coords$y
  map <- map %>% relocate(y) %>% relocate(x)
  res[[i]] <- map
}

rctd_plot <- ggplot() +
  geom_point(
    data = res[[1]], mapping = aes(x = x, y = y, alpha = `T cells`),
    color = "green", size = 0.5
  ) +
  geom_point(
    data = res[[1]], mapping = aes(x = x, y = y, alpha = `B cells`),
    color = "red", size = 0.5
  ) +
  geom_point(
    data = res[[1]], mapping = aes(x = x, y = y, alpha = `Epithelial Cells`),
    color = "cyan", size = 0.5
  ) +
  theme_minimal() +
  scale_alpha_continuous(range = c(0, 1), breaks = c(0, 0.1, 1)) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    axis.line = element_blank(),
    text = element_blank(),
    plot.title = element_text(size = 30),
    plot.title.position = "panel",
    axis.title = element_blank()
  ) +
  coord_fixed()

print(rctd_plot)

output_file <- Sys.getenv("FLU_RCTD_VIS_OUTPUT", unset = "")
if (nzchar(output_file)) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_file, plot = rctd_plot)
  message("Saved RCTD visualization to: ", output_file)
}
