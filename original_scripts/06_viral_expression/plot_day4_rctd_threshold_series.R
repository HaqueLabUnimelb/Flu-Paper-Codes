#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

project_dir <- "/data/projects/punim2394/lung_curio_mouse_pr8"
analysis_dir <- file.path(project_dir, "viralgeneassess")
output_dir <- file.path(analysis_dir, "thresholded_celltype_barline")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_umi_thresholds <- c(seq(1L, 5L, by = 1L), seq(10L, 50L, by = 5L))
samples <- data.frame(
  sample_id = c("2_2", "2_3", "4_2", "4_3"),
  day = c("Day 4", "Day 4", "Day 10", "Day 10"),
  input_file = file.path(
    project_dir, "processed_data",
    c("annotated_lungST_2_2_QCed.rds", "annotated_lungST_2_3_QCed.rds",
      "annotated_lungST_4_2_QCed.rds", "annotated_lungST_4_3_QCed.rds")
  ),
  stringsAsFactors = FALSE
)
stopifnot(all(file.exists(samples$input_file)))

get_counts <- function(object) {
  assay_name <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
  count_layers <- grep("^counts($|\\.)", Layers(object[[assay_name]]), value = TRUE)
  if (length(count_layers) != 1L) {
    stop("Expected exactly one RNA counts layer, found: ", paste(count_layers, collapse = ", "))
  }
  LayerData(object, assay = assay_name, layer = count_layers)
}

detect_pr8_features <- function(feature_names) {
  tagged <- grepl(
    "(^|[-_.:|/])pr8([-_.:|/]|$)|puerto[ _/-]*rico|influenza[ _-]*a|h1n1",
    feature_names, ignore.case = TRUE, perl = TRUE
  ) | grepl(
    "(^|[-_.:|])iav([-_.:|]|$)",
    feature_names, ignore.case = TRUE, perl = TRUE
  )
  canonical <- toupper(feature_names) %in% c(
    "PB2", "PB1", "PA", "HA", "NP", "NA",
    "M", "M1", "M2", "NS", "NS1", "NS2", "NEP"
  )
  sort(unique(feature_names[tagged | canonical]))
}

get_spatial_coordinates <- function(object) {
  reductions <- Reductions(object)
  spatial_reduction <- reductions[tolower(reductions) == "spatial"]
  if (length(spatial_reduction) == 0L) {
    spatial_reduction <- grep("spatial", reductions, ignore.case = TRUE, value = TRUE)
  }
  if (length(spatial_reduction) == 0L) {
    stop("No spatial reduction found. Available reductions: ", paste(reductions, collapse = ", "))
  }
  coordinates <- Embeddings(object, reduction = spatial_reduction[[1L]])
  if (ncol(coordinates) < 2L) {
    stop("Spatial reduction has fewer than two dimensions.")
  }
  coordinates[, 1:2, drop = FALSE]
}

accepted_classes <- c("singlet", "doublet_uncertain", "doublet_certain")
spot_tables <- list()
spatial_tables <- list()

for (index in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[[index]]
  message("Reading sample ", sample_id, ": ", samples$input_file[[index]])
  object <- readRDS(samples$input_file[[index]])

  counts <- get_counts(object)
  pr8_features <- detect_pr8_features(rownames(counts))
  if (length(pr8_features) == 0L) {
    stop("No PR8 features detected in sample ", sample_id, ".")
  }
  pr8_total_umi <- Matrix::colSums(counts[pr8_features, , drop = FALSE])
  names(pr8_total_umi) <- colnames(counts)
  coordinates <- get_spatial_coordinates(object)
  if (!all(names(pr8_total_umi) %in% rownames(coordinates))) {
    stop("Counts contain spots without spatial coordinates in sample ", sample_id, ".")
  }
  coordinates <- coordinates[names(pr8_total_umi), , drop = FALSE]
  spatial_tables[[sample_id]] <- data.frame(
    sample_id = sample_id,
    day = samples$day[[index]],
    spot = names(pr8_total_umi),
    spatial_1 = coordinates[, 1L],
    spatial_2 = coordinates[, 2L],
    pr8_total_raw_umi = as.numeric(pr8_total_umi),
    stringsAsFactors = FALSE
  )

  if (samples$day[[index]] == "Day 4") {
    metadata <- object[[]]
    required_metadata <- c("spot_class", "first_type")
    missing_metadata <- setdiff(required_metadata, colnames(metadata))
    if (length(missing_metadata) > 0L) {
      stop("Missing RCTD metadata in ", sample_id, ": ", paste(missing_metadata, collapse = ", "))
    }
    if (!setequal(rownames(metadata), names(pr8_total_umi))) {
      stop("Metadata and counts have different spot names in sample ", sample_id, ".")
    }
    metadata <- metadata[names(pr8_total_umi), , drop = FALSE]
    accepted <- metadata$spot_class %in% accepted_classes &
      !is.na(metadata$first_type) & nzchar(as.character(metadata$first_type))
    spot_tables[[sample_id]] <- data.frame(
      sample_id = sample_id,
      cell_type = as.character(metadata$first_type[accepted]),
      pr8_total_raw_umi = as.numeric(pr8_total_umi[accepted]),
      stringsAsFactors = FALSE
    )
  }
  rm(object, counts)
  invisible(gc())
}

day4_spots <- do.call(rbind, spot_tables)
all_spatial <- do.call(rbind, spatial_tables)
cell_totals <- aggregate(
  rep(1L, nrow(day4_spots)),
  by = list(cell_type = day4_spots$cell_type),
  FUN = sum
)
colnames(cell_totals)[[2L]] <- "total_spots"
cell_totals <- cell_totals[order(cell_totals$total_spots, decreasing = TRUE), ]
cell_type_levels <- cell_totals$cell_type

summaries <- lapply(raw_umi_thresholds, function(raw_umi_threshold) {
  infected <- day4_spots$pr8_total_raw_umi >= raw_umi_threshold
  result <- aggregate(
    cbind(total_spots = rep(1L, nrow(day4_spots)), infected_spots = as.integer(infected)) ~ cell_type,
    data = day4_spots,
    FUN = sum
  )
  result$uninfected_spots <- result$total_spots - result$infected_spots
  result$infected_percent <- 100 * result$infected_spots / result$total_spots
  result$raw_umi_threshold <- raw_umi_threshold
  result <- result[match(cell_type_levels, result$cell_type), ]
  result
})
names(summaries) <- as.character(raw_umi_thresholds)
combined_summary <- do.call(rbind, summaries)

write.table(
  combined_summary,
  file.path(output_dir, "Day4_RCTD_celltype_raw_UMI_threshold_series.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

common_bar_max <- 1.06 * max(combined_summary$infected_spots)
positive_spatial_umi <- all_spatial$pr8_total_raw_umi[all_spatial$pr8_total_raw_umi > 0]
spatial_colour_cap_raw <- max(
  as.numeric(quantile(positive_spatial_umi, 0.99, names = FALSE)),
  max(raw_umi_thresholds) + 1
)
spatial_colour_cap_log1p <- log1p(spatial_colour_cap_raw)

for (threshold_index in seq_along(raw_umi_thresholds)) {
  raw_umi_threshold <- raw_umi_thresholds[[threshold_index]]
  threshold_label <- sprintf("%02d", raw_umi_threshold)
  threshold_dir <- file.path(output_dir, paste0("raw_UMI_threshold_", threshold_label))
  dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)
  plot_data <- summaries[[threshold_index]]
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_levels)

  percentage_plot <- ggplot(
    plot_data, aes(x = cell_type, y = infected_percent, group = 1)
  ) +
    geom_line(colour = "#7F0000", linewidth = 0.9) +
    geom_point(
      shape = 21, size = 3.0, stroke = 0.7,
      colour = "#7F0000", fill = "#EF3B2C"
    ) +
    geom_text(
      aes(label = sprintf("%.1f%%", infected_percent)),
      vjust = -0.75, size = 3.2, colour = "#B2182B"
    ) +
    scale_y_continuous(
      limits = c(0, 100), breaks = seq(0, 100, by = 20),
      labels = function(value) paste0(value, "%"),
      expand = expansion(mult = c(0.01, 0.06))
    ) +
    labs(x = NULL, y = "Threshold-positive spots (%)") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.line.x = element_blank(), plot.margin = margin(8, 15, 0, 12)
    )

  bar_plot <- ggplot(plot_data, aes(x = cell_type, y = infected_spots)) +
    geom_col(width = 0.68, fill = "#D73027", colour = "grey25", linewidth = 0.3) +
    scale_y_continuous(
      limits = c(0, common_bar_max), labels = scales::comma,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = "Number of threshold-positive spots") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = margin(0, 15, 8, 12)
    )

  combined_plot <- percentage_plot / bar_plot +
    plot_layout(heights = c(1.25, 3.0)) +
    plot_annotation(
      title = paste0("Day 4 PR8 signal by RCTD cell type — raw UMI threshold ", raw_umi_threshold),
      subtitle = paste0(
        "Replicates 2_2 and 2_3 pooled | accepted RCTD spots only | positive = ",
        "total PR8 raw UMI >= ", raw_umi_threshold
      ),
      caption = paste0(
        "Bars show threshold-positive spot counts. The aligned line and labels show ",
        "the threshold-positive percentage within each cell type."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11),
        plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
        plot.margin = margin(12, 12, 12, 12)
      )
    )

  output_file <- file.path(
    threshold_dir,
    paste0("Day4_RCTD_celltype_barline_raw_UMI_threshold_", threshold_label, ".png")
  )
  ggsave(
    output_file, combined_plot,
    width = 13, height = 9, units = "in", dpi = 300, bg = "white"
  )
  message("Saved: ", output_file)

  write.table(
    plot_data,
    file.path(threshold_dir, "Day4_RCTD_celltype_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  for (sample_index in seq_len(nrow(samples))) {
    sample_id <- samples$sample_id[[sample_index]]
    spatial_data <- spatial_tables[[sample_id]]
    threshold_positive <- spatial_data$pr8_total_raw_umi >= raw_umi_threshold
    positive_data <- spatial_data[threshold_positive, , drop = FALSE]
    positive_data <- positive_data[order(positive_data$pr8_total_raw_umi), , drop = FALSE]
    positive_data$log1p_pr8_total_umi <- log1p(positive_data$pr8_total_raw_umi)
    n_positive <- nrow(positive_data)

    spatial_plot <- ggplot(
      spatial_data, aes(x = spatial_1, y = spatial_2)
    ) +
      geom_point(size = 0.077, colour = "grey90")
    if (n_positive > 0L) {
      spatial_plot <- spatial_plot +
        geom_point(
          data = positive_data,
          aes(colour = log1p_pr8_total_umi),
          size = 0.091, alpha = 1
        )
    }
    spatial_plot <- spatial_plot +
      coord_fixed() +
      scale_colour_viridis_c(
        option = "magma",
        direction = -1,
        limits = c(log1p(raw_umi_threshold), spatial_colour_cap_log1p),
        oob = scales::squish,
        labels = scales::label_number(accuracy = 0.1),
        name = "log1p(total\nPR8 UMI)"
      ) +
      theme_void(base_size = 11) +
      theme(
        legend.position = "right",
        plot.margin = margin(4, 4, 4, 4)
      )

    spatial_output <- file.path(
      threshold_dir,
      sprintf("PR8_viral_gene_expression_%s.png", sample_id)
    )
    ggsave(
      spatial_output, spatial_plot,
      width = 8, height = 7, units = "in", dpi = 300, bg = "white"
    )
    message("Saved: ", spatial_output)
  }
}

legacy_files <- c(
  "Day4_RCTD_celltype_barline_log1p_cutoff_0p2.png",
  "Day4_RCTD_celltype_barline_log1p_cutoff_0p3.png",
  "Day4_RCTD_celltype_barline_log1p_cutoff_0p5.png",
  "Day4_RCTD_celltype_barline_log1p_cutoff_0p8.png",
  "Day4_RCTD_celltype_barline_log1p_cutoff_0p9.png",
  "Day4_RCTD_celltype_barline_log1p_cutoff_1p0.png",
  "cutoff_to_raw_umi_mapping.tsv",
  "Day4_RCTD_celltype_threshold_series.tsv"
)
legacy_paths <- file.path(output_dir, legacy_files)
legacy_existing <- legacy_paths[file.exists(legacy_paths)]
if (length(legacy_existing) > 0L) {
  unlink(legacy_existing)
  message("Removed superseded log1p-cutoff outputs: ", length(legacy_existing))
}

writeLines(
  c(
    "PR8 raw-UMI threshold series",
    "",
    "Thresholds: total PR8 raw UMI per spot >= 1, 2, 3, 4, 5, 10, 15, ..., 50.",
    "Each raw_UMI_threshold_XX subfolder contains one Day 4 RCTD cell-type bar/line plot, its TSV summary, and four spatial plots.",
    "The RCTD plot pools accepted spots from samples 2_2 and 2_3.",
    "Spatial plots include samples 2_2, 2_3, 4_2, and 4_3; grey spots are below threshold and coloured spots meet the raw-UMI threshold.",
    "The cell-type order, bar y-axis limits, and spatial colour upper cap are fixed across thresholds for direct comparison.",
    sprintf("Shared spatial colour upper cap: log1p(%.0f total PR8 raw UMI) = %.3f.", spatial_colour_cap_raw, spatial_colour_cap_log1p),
    "Spatial colour is log1p(total PR8 UMI), the colour scale is reversed, point sizes are halved, and only legend text is retained."
  ),
  file.path(output_dir, "README.txt")
)
