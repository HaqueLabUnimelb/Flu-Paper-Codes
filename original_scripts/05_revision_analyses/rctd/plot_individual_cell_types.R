library(spacexr)
library(ggplot2)

input_file <- paste0(
  "/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/",
  "RCTD/Day4mouseRCTDoutmin30_multi.rds"
)
output_dir <- paste0(
  "/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/",
  "RCTD/Day4mouseRCTDoutmin30_multi_cell_type_plots"
)
replicate_names <- c("puck1_lungST_2_2", "puck2_lungST_2_3")

safe_filename <- function(x) {
  gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

rctd <- readRDS(input_file)
if (!inherits(rctd, "RCTD.replicates")) {
  stop("Expected an RCTD.replicates object, got: ", paste(class(rctd), collapse = ", "))
}
if (length(rctd@RCTD.reps) != length(replicate_names)) {
  stop("Unexpected number of RCTD replicates: ", length(rctd@RCTD.reps))
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
