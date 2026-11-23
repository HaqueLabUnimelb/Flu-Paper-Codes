#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

project_dir <- "/data/projects/punim2394/lung_curio_mouse_pr8"
output_dir <- file.path(project_dir, "viralgeneassess")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- data.frame(
  sample_id = c("2_2", "2_3", "4_2", "4_3"),
  day = c("Day 4", "Day 4", "Day 10", "Day 10"),
  input_file = file.path(
    project_dir,
    "processed_data",
    c(
      "annotated_lungST_2_2_QCed.rds",
      "annotated_lungST_2_3_QCed.rds",
      "annotated_lungST_4_2_QCed.rds",
      "annotated_lungST_4_3_QCed.rds"
    )
  ),
  stringsAsFactors = FALSE
)

stopifnot(all(file.exists(samples$input_file)))

get_counts <- function(object) {
  assay_name <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
  assay_object <- object[[assay_name]]
  layer_names <- Layers(assay_object)
  count_layers <- grep("^counts($|\\.)", layer_names, value = TRUE)

  if (length(count_layers) != 1L) {
    stop(
      "Expected exactly one counts layer in assay '", assay_name,
      "', but found: ", paste(count_layers, collapse = ", ")
    )
  }

  counts <- LayerData(object, assay = assay_name, layer = count_layers)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("The counts matrix must have feature and spot names.")
  }

  list(matrix = counts, assay = assay_name, layer = count_layers)
}

detect_pr8_features <- function(feature_names) {
  tagged <- grepl(
    "(^|[-_.:|/])pr8([-_.:|/]|$)|puerto[ _/-]*rico|influenza[ _-]*a|h1n1",
    feature_names,
    ignore.case = TRUE,
    perl = TRUE
  ) | grepl(
    "(^|[-_.:|])iav([-_.:|]|$)",
    feature_names,
    ignore.case = TRUE,
    perl = TRUE
  )

  canonical_products <- c(
    "PB2", "PB1", "PA", "HA", "NP", "NA",
    "M", "M1", "M2", "NS", "NS1", "NS2", "NEP"
  )
  canonical <- toupper(feature_names) %in% canonical_products

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

  list(matrix = coordinates[, 1:2, drop = FALSE], reduction = spatial_reduction[[1L]])
}

plot_data <- list()
genes_by_sample <- list()
gene_summaries <- list()
sample_summaries <- list()
processing_notes <- character()

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[[i]]
  message("Reading sample ", sample_id, ": ", samples$input_file[[i]])
  object <- readRDS(samples$input_file[[i]])

  counts_info <- get_counts(object)
  counts <- counts_info$matrix
  pr8_features <- detect_pr8_features(rownames(counts))
  if (length(pr8_features) == 0L) {
    stop(
      "No PR8 features detected in sample ", sample_id,
      ". Inspect feature names before changing the detection rules."
    )
  }

  total_pr8_umi <- Matrix::colSums(counts[pr8_features, , drop = FALSE])
  names(total_pr8_umi) <- colnames(counts)

  coordinate_info <- get_spatial_coordinates(object)
  coordinates <- coordinate_info$matrix
  missing_coordinates <- setdiff(names(total_pr8_umi), rownames(coordinates))
  if (length(missing_coordinates) > 0L) {
    stop(
      "Sample ", sample_id, " has ", length(missing_coordinates),
      " spots without spatial coordinates."
    )
  }
  coordinates <- coordinates[names(total_pr8_umi), , drop = FALSE]

  sample_plot_data <- data.frame(
    sample_id = sample_id,
    day = samples$day[[i]],
    spot = rownames(coordinates),
    spatial_1 = coordinates[, 1L],
    spatial_2 = coordinates[, 2L],
    pr8_total_umi = as.numeric(total_pr8_umi),
    stringsAsFactors = FALSE
  )
  sample_plot_data$pr8_viral_gene_expression <- log1p(sample_plot_data$pr8_total_umi)
  plot_data[[sample_id]] <- sample_plot_data

  feature_matrix <- counts[pr8_features, , drop = FALSE]
  genes_by_sample[[sample_id]] <- data.frame(
    sample_id = sample_id,
    day = samples$day[[i]],
    feature = pr8_features,
    assay = counts_info$assay,
    layer = counts_info$layer,
    stringsAsFactors = FALSE
  )
  gene_summaries[[sample_id]] <- data.frame(
    sample_id = sample_id,
    day = samples$day[[i]],
    feature = pr8_features,
    total_umi = as.numeric(Matrix::rowSums(feature_matrix)),
    spots_detected = as.integer(Matrix::rowSums(feature_matrix > 0)),
    percent_spots_detected = 100 * as.numeric(Matrix::rowSums(feature_matrix > 0)) / ncol(feature_matrix),
    stringsAsFactors = FALSE
  )
  sample_summaries[[sample_id]] <- data.frame(
    sample_id = sample_id,
    day = samples$day[[i]],
    spots = length(total_pr8_umi),
    pr8_features = length(pr8_features),
    total_pr8_umi = sum(total_pr8_umi),
    spots_pr8_positive = sum(total_pr8_umi > 0),
    percent_spots_pr8_positive = 100 * mean(total_pr8_umi > 0),
    mean_pr8_umi_per_spot = mean(total_pr8_umi),
    median_pr8_umi_per_spot = median(total_pr8_umi),
    p90_pr8_umi_per_spot = as.numeric(quantile(total_pr8_umi, 0.90, names = FALSE)),
    p99_pr8_umi_per_spot = as.numeric(quantile(total_pr8_umi, 0.99, names = FALSE)),
    max_pr8_umi_per_spot = max(total_pr8_umi),
    stringsAsFactors = FALSE
  )
  processing_notes <- c(
    processing_notes,
    sprintf(
      "%s (%s): assay=%s, layer=%s, reduction=%s, PR8 features=%s",
      sample_id,
      samples$day[[i]],
      counts_info$assay,
      counts_info$layer,
      coordinate_info$reduction,
      paste(pr8_features, collapse = ", ")
    )
  )

  rm(object, counts, feature_matrix)
  invisible(gc())
}

genes_by_sample <- do.call(rbind, genes_by_sample)
gene_summaries <- do.call(rbind, gene_summaries)
sample_summaries <- do.call(rbind, sample_summaries)
all_plot_data <- do.call(rbind, plot_data)
pr8_feature_union <- sort(unique(genes_by_sample$feature))

write.table(
  genes_by_sample,
  file.path(output_dir, "pr8_features_by_sample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(pr8_feature_union, file.path(output_dir, "pr8_feature_union.txt"))
write.table(
  gene_summaries,
  file.path(output_dir, "pr8_feature_expression_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  sample_summaries,
  file.path(output_dir, "pr8_sample_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

positive_values <- all_plot_data$pr8_total_umi[all_plot_data$pr8_total_umi > 0]
common_cap_raw <- if (length(positive_values) == 0L) {
  1
} else {
  as.numeric(quantile(positive_values, 0.99, names = FALSE))
}
common_cap_log1p <- log1p(max(common_cap_raw, 1))

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[[i]]
  data_to_plot <- plot_data[[sample_id]]
  data_to_plot <- data_to_plot[order(data_to_plot$pr8_viral_gene_expression), ]
  summary_row <- sample_summaries[sample_summaries$sample_id == sample_id, ]

  plot_object <- ggplot(
    data_to_plot,
    aes(x = spatial_1, y = spatial_2, colour = pr8_viral_gene_expression)
  ) +
    geom_point(size = 0.55, alpha = 1) +
    coord_fixed() +
    scale_colour_viridis_c(
      option = "magma",
      limits = c(0, common_cap_log1p),
      oob = scales::squish,
      name = "PR8 viral gene expression\nlog1p(total UMI)"
    ) +
    labs(
      title = paste0("PR8 viral gene expression — sample ", sample_id),
      subtitle = sprintf(
        "%s | %d PR8 features | PR8-positive spots: %d/%d (%.2f%%)",
        samples$day[[i]],
        summary_row$pr8_features,
        summary_row$spots_pr8_positive,
        summary_row$spots,
        summary_row$percent_spots_pr8_positive
      ),
      caption = sprintf(
        "Raw RNA counts summed per spot; colour shown as log1p and capped at the shared positive-spot 99th percentile (%.0f UMI).",
        common_cap_raw
      )
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8, colour = "grey35", hjust = 0),
      legend.position = "right",
      plot.margin = margin(12, 12, 12, 12)
    )

  output_png <- file.path(
    output_dir,
    sprintf("PR8_viral_gene_expression_%s.png", sample_id)
  )
  ggsave(output_png, plot_object, width = 8, height = 7, units = "in", dpi = 300, bg = "white")
  message("Saved: ", output_png)
}

report_lines <- c(
  "PR8 viral gene expression assessment",
  "",
  "Definition: for each spatial spot, raw RNA-assay UMI counts are summed across all detected PR8 viral features. The plotted value is log1p(total PR8 UMI).",
  sprintf("Shared plot cap: positive-spot 99th percentile = %.0f raw UMI; log1p = %.6f.", common_cap_raw, common_cap_log1p),
  "",
  "PR8 feature union:",
  paste(pr8_feature_union, collapse = ", "),
  "",
  "Per-sample processing:",
  processing_notes
)
writeLines(report_lines, file.path(output_dir, "README_results.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Finished PR8 viral gene expression assessment.")
