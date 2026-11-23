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
  sample_id = c("2_2", "2_3"),
  input_file = file.path(
    project_dir,
    "processed_data",
    c(
      "annotated_lungST_2_2_QCed.rds",
      "annotated_lungST_2_3_QCed.rds"
    )
  ),
  stringsAsFactors = FALSE
)

stopifnot(all(file.exists(samples$input_file)))

get_counts <- function(object) {
  assay_name <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
  layer_names <- Layers(object[[assay_name]])
  count_layers <- grep("^counts($|\\.)", layer_names, value = TRUE)
  if (length(count_layers) != 1L) {
    stop(
      "Expected exactly one counts layer in assay '", assay_name,
      "', but found: ", paste(count_layers, collapse = ", ")
    )
  }
  counts <- LayerData(object, assay = assay_name, layer = count_layers)
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

accepted_classes <- c("singlet", "doublet_uncertain", "doublet_certain")
spot_tables <- list()
processing_notes <- character()

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[[i]]
  message("Reading Day 4 sample ", sample_id, ": ", samples$input_file[[i]])
  object <- readRDS(samples$input_file[[i]])
  metadata <- object[[]]

  required_columns <- c("spot_class", "first_type")
  missing_columns <- setdiff(required_columns, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      "Sample ", sample_id, " is missing RCTD metadata columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  counts_info <- get_counts(object)
  counts <- counts_info$matrix
  pr8_features <- detect_pr8_features(rownames(counts))
  if (length(pr8_features) == 0L) {
    stop("No PR8 features detected in sample ", sample_id, ".")
  }

  total_pr8_umi <- Matrix::colSums(counts[pr8_features, , drop = FALSE])
  names(total_pr8_umi) <- colnames(counts)
  if (!setequal(rownames(metadata), names(total_pr8_umi))) {
    stop("Metadata and RNA counts contain different spot names in sample ", sample_id, ".")
  }
  metadata <- metadata[names(total_pr8_umi), , drop = FALSE]

  all_spots <- data.frame(
    sample_id = sample_id,
    spot = rownames(metadata),
    spot_class = as.character(metadata$spot_class),
    cell_type = as.character(metadata$first_type),
    pr8_total_umi = as.numeric(total_pr8_umi),
    infected = as.numeric(total_pr8_umi) > 0,
    stringsAsFactors = FALSE
  )
  accepted <- all_spots$spot_class %in% accepted_classes &
    !is.na(all_spots$cell_type) & nzchar(all_spots$cell_type)
  annotated_spots <- all_spots[accepted, , drop = FALSE]
  if (nrow(annotated_spots) == 0L) {
    stop("No accepted RCTD-annotated spots found in sample ", sample_id, ".")
  }
  spot_tables[[sample_id]] <- annotated_spots

  processing_notes <- c(
    processing_notes,
    sprintf(
      paste0(
        "%s: total spots=%d; accepted RCTD spots=%d; rejected/excluded=%d; ",
        "PR8-positive accepted spots=%d (%.2f%%); PR8 features=%s"
      ),
      sample_id,
      nrow(all_spots),
      nrow(annotated_spots),
      nrow(all_spots) - nrow(annotated_spots),
      sum(annotated_spots$infected),
      100 * mean(annotated_spots$infected),
      paste(pr8_features, collapse = ", ")
    )
  )

  rm(object, metadata, counts, all_spots, annotated_spots)
  invisible(gc())
}

day4_spots <- do.call(rbind, spot_tables)

summarise_infection <- function(data, grouping_variables) {
  formula_text <- paste(
    "cbind(total_spots = rep(1L, nrow(data)), infected_spots = as.integer(infected)) ~",
    paste(grouping_variables, collapse = " + ")
  )
  summary_table <- aggregate(
    as.formula(formula_text),
    data = data,
    FUN = sum
  )
  summary_table$uninfected_spots <- summary_table$total_spots - summary_table$infected_spots
  summary_table$infected_percent <- 100 * summary_table$infected_spots / summary_table$total_spots
  summary_table
}

pooled_summary <- summarise_infection(day4_spots, "cell_type")
pooled_summary <- pooled_summary[order(pooled_summary$total_spots, decreasing = TRUE), ]
replicate_summary <- summarise_infection(day4_spots, c("sample_id", "cell_type"))
replicate_summary <- replicate_summary[order(replicate_summary$sample_id, -replicate_summary$total_spots), ]

stopifnot(
  sum(pooled_summary$total_spots) == nrow(day4_spots),
  all(pooled_summary$infected_spots <= pooled_summary$total_spots),
  sum(pooled_summary$infected_spots) == sum(day4_spots$infected)
)

write.table(
  pooled_summary,
  file.path(output_dir, "Day4_RCTD_celltype_infection_pooled.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  replicate_summary,
  file.path(output_dir, "Day4_RCTD_celltype_infection_by_replicate.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cell_type_levels <- pooled_summary$cell_type
plot_data <- rbind(
  data.frame(
    cell_type = pooled_summary$cell_type,
    bar_type = "All annotated spots",
    number_of_spots = pooled_summary$total_spots,
    infected_percent = pooled_summary$infected_percent,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cell_type = pooled_summary$cell_type,
    bar_type = "Infected spots",
    number_of_spots = pooled_summary$infected_spots,
    infected_percent = pooled_summary$infected_percent,
    stringsAsFactors = FALSE
  )
)
plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_levels)
plot_data$bar_type <- factor(
  plot_data$bar_type,
  levels = c("All annotated spots", "Infected spots")
)
plot_data$label <- ifelse(
  plot_data$bar_type == "Infected spots",
  sprintf("%.1f%%", plot_data$infected_percent),
  ""
)

bar_plot <- ggplot(
  plot_data,
  aes(x = cell_type, y = number_of_spots, fill = bar_type)
) +
  geom_col(
    position = position_dodge(width = 0.80),
    width = 0.72,
    colour = "grey25",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.80),
    vjust = -0.35,
    size = 3.0,
    colour = "#B2182B"
  ) +
  scale_fill_manual(
    values = c(
      "All annotated spots" = "#BDBDBD",
      "Infected spots" = "#D73027"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.13))
  ) +
  labs(
    title = "Day 4 PR8 infection by RCTD cell type",
    subtitle = paste0(
      "Replicates 2_2 and 2_3 pooled | accepted RCTD spots only | ",
      "infected = total PR8 UMI > 0"
    ),
    x = NULL,
    y = "Number of spots",
    caption = paste0(
      "Cell type is the dominant RCTD label (first_type). ",
      "Red labels show the infected percentage within each cell type."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(12, 16, 12, 12)
  )

plot_file <- file.path(output_dir, "Day4_RCTD_celltype_infected_spots_barplot.png")
ggsave(plot_file, bar_plot, width = 13, height = 8, units = "in", dpi = 300, bg = "white")

report_lines <- c(
  "Day 4 RCTD cell-type PR8 infection assessment",
  "",
  "Replicates: 2_2 and 2_3, pooled at the spot level.",
  paste0("Included RCTD classes: ", paste(accepted_classes, collapse = ", "), "."),
  "Excluded RCTD class: reject.",
  "Cell-type assignment: first_type (dominant RCTD cell type).",
  "Infected definition: summed raw UMI across detected PR8 viral features > 0.",
  sprintf(
    "Pooled accepted spots: %d; infected: %d (%.2f%%).",
    nrow(day4_spots),
    sum(day4_spots$infected),
    100 * mean(day4_spots$infected)
  ),
  "",
  "Per-replicate processing:",
  processing_notes
)
writeLines(report_lines, file.path(output_dir, "README_Day4_RCTD_celltype_infection.txt"))

message("Saved plot: ", plot_file)
message("Pooled accepted spots: ", nrow(day4_spots))
message("Pooled infected accepted spots: ", sum(day4_spots$infected))
