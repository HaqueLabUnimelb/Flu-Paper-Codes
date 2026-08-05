#!/usr/bin/env Rscript
# Purpose: render one spatial proportion plot per cell type and day-4 replicate.
# Input: FLU_RCTD_DAY4_FULL_RDS.
# Output directory: FLU_RCTD_CELL_PLOT_OUTPUT_DIR.

suppressPackageStartupMessages({
  library(spacexr)
  library(ggplot2)
})

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
output_dir <- env_path(
  "FLU_RCTD_CELL_PLOT_OUTPUT_DIR",
  file.path(repo_dir, "outputs", "revision", "rctd", "day4_cell_type_plots")
)
replicate_names <- trimws(strsplit(
  Sys.getenv("FLU_RCTD_DAY4_REPLICATE_NAMES", "puck1_lungST_2_2,puck2_lungST_2_3"),
  ",",
  fixed = TRUE
)[[1L]])

safe_filename <- function(x) {
  gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

rctd <- readRDS(input_file)
if (!inherits(rctd, "RCTD.replicates")) {
  stop("Expected an RCTD.replicates object, got: ", paste(class(rctd), collapse = ", "))
}
if (length(replicate_names) != length(rctd@RCTD.reps)) {
  stop("FLU_RCTD_DAY4_REPLICATE_NAMES does not match the replicate count")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summary_rows <- list()

for (i in seq_along(rctd@RCTD.reps)) {
  replicate <- rctd@RCTD.reps[[i]]
  weights <- as.matrix(spacexr::normalize_weights(replicate@results$weights))
  coords <- as.data.frame(replicate@spatialRNA@coords)

  if (!identical(rownames(weights), rownames(coords))) {
    missing_coords <- setdiff(rownames(weights), rownames(coords))
    if (length(missing_coords) > 0L) {
      stop("Missing coordinates for ", length(missing_coords), " spots in replicate ", i)
    }
    coords <- coords[rownames(weights), , drop = FALSE]
  }
  colnames(coords)[1:2] <- c("x", "y")

  replicate_dir <- file.path(output_dir, replicate_names[[i]])
  dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)

  for (cell_type in colnames(weights)) {
    plot_data <- data.frame(
      x = coords$x,
      y = coords$y,
      weight = weights[, cell_type],
      stringsAsFactors = FALSE
    )
    plot_data <- plot_data[order(plot_data$weight), , drop = FALSE]

    plot_object <- ggplot(plot_data, aes(x = x, y = y, colour = weight)) +
      geom_point(size = 0.35, alpha = 1) +
      scale_colour_viridis_c(
        option = "magma",
        limits = c(0, 1),
        oob = scales::squish,
        name = "RCTD\nproportion"
      ) +
      coord_fixed() +
      labs(title = paste(replicate_names[[i]], cell_type, sep = " — ")) +
      theme_void(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA)
      )

    output_file <- file.path(
      replicate_dir,
      paste0(sprintf("%02d", match(cell_type, colnames(weights))), "_", safe_filename(cell_type), ".png")
    )
    ggsave(
      filename = output_file,
      plot = plot_object,
      width = 7,
      height = 6,
      units = "in",
      dpi = 300,
      bg = "white"
    )

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      replicate = replicate_names[[i]],
      cell_type = cell_type,
      spots = nrow(plot_data),
      min_weight = min(plot_data$weight),
      mean_weight = mean(plot_data$weight),
      max_weight = max(plot_data$weight),
      output_file = output_file,
      stringsAsFactors = FALSE
    )
    message("Saved: ", output_file)
  }
}

summary_table <- do.call(rbind, summary_rows)
write.csv(
  summary_table,
  file.path(output_dir, "cell_type_plot_summary.csv"),
  row.names = FALSE
)
message("Completed ", nrow(summary_table), " cell-type plots in: ", output_dir)
