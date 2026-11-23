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

segment_levels <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

get_counts <- function(object) {
  assay_name <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
  count_layers <- grep(
    "^counts($|\\.)", Layers(object[[assay_name]]), value = TRUE
  )
  if (length(count_layers) != 1L) {
    stop(
      "Expected exactly one counts layer in assay '", assay_name,
      "', but found: ", paste(count_layers, collapse = ", ")
    )
  }
  LayerData(object, assay = assay_name, layer = count_layers)
}

map_segment <- function(feature_names) {
  products <- toupper(feature_names)
  segment <- rep(NA_character_, length(products))
  canonical <- c("PB2", "PB1", "PA", "HA", "NP", "NA")
  segment[products %in% canonical] <- products[products %in% canonical]
  segment[products %in% c("M", "M1", "M2")] <- "M"
  segment[products %in% c("NS", "NS1", "NS2", "NEP")] <- "NS"
  segment
}

spot_tables <- list()
mapping_tables <- list()
summary_tables <- list()

for (index in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[[index]]
  message("Reading sample ", sample_id, ": ", samples$input_file[[index]])
  object <- readRDS(samples$input_file[[index]])
  counts <- get_counts(object)
  segment_assignment <- map_segment(rownames(counts))
  viral_features <- rownames(counts)[!is.na(segment_assignment)]
  viral_segments <- segment_assignment[!is.na(segment_assignment)]
  if (length(viral_features) == 0L) {
    stop("No canonical PR8 segment features found in sample ", sample_id, ".")
  }

  mapping_tables[[sample_id]] <- data.frame(
    sample_id = sample_id,
    day = samples$day[[index]],
    segment = segment_levels,
    input_features = vapply(
      segment_levels,
      function(segment_name) {
        paste(viral_features[viral_segments == segment_name], collapse = "+")
      },
      character(1L)
    ),
    stringsAsFactors = FALSE
  )
  mapping_tables[[sample_id]]$segment_present <-
    nzchar(mapping_tables[[sample_id]]$input_features)

  total_viral_umi <- as.numeric(
    Matrix::colSums(counts[viral_features, , drop = FALSE])
  )
  spot_tables[[sample_id]] <- data.frame(
    sample_id = sample_id,
    day = samples$day[[index]],
    total_pr8_8segment_umi = total_viral_umi,
    log1p_total_pr8_8segment_umi = log1p(total_viral_umi),
    stringsAsFactors = FALSE
  )
  summary_tables[[sample_id]] <- data.frame(
    sample_id = sample_id,
    day = samples$day[[index]],
    spots = length(total_viral_umi),
    segments_present = length(unique(viral_segments)),
    zero_spots = sum(total_viral_umi == 0),
    positive_spots = sum(total_viral_umi > 0),
    positive_percent = 100 * mean(total_viral_umi > 0),
    median_raw_umi = median(total_viral_umi),
    p90_raw_umi = as.numeric(quantile(total_viral_umi, 0.90, names = FALSE)),
    p99_raw_umi = as.numeric(quantile(total_viral_umi, 0.99, names = FALSE)),
    max_raw_umi = max(total_viral_umi),
    stringsAsFactors = FALSE
  )

  rm(object, counts)
  invisible(gc())
}

plot_data <- do.call(rbind, spot_tables)
mapping_table <- do.call(rbind, mapping_tables)
summary_table <- do.call(rbind, summary_tables)
plot_data$sample_id <- factor(plot_data$sample_id, levels = samples$sample_id)

write.table(
  mapping_table,
  file.path(output_dir, "PR8_8segment_feature_mapping_by_sample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  summary_table,
  file.path(output_dir, "PR8_8segment_spot_distribution_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# Count spots at each exact integer value of total PR8 raw UMI. This table and
# plot pool the four samples; sample-level values remain available in plot_data.
bar_data <- as.data.frame(table(plot_data$total_pr8_8segment_umi))
names(bar_data) <- c("total_pr8_8segment_umi", "number_of_spots")
bar_data$total_pr8_8segment_umi <- as.numeric(
  as.character(bar_data$total_pr8_8segment_umi)
)
bar_data$number_of_spots <- as.integer(bar_data$number_of_spots)

write.table(
  bar_data,
  file.path(output_dir, "PR8_8segment_total_UMI_spot_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

bar_plot <- ggplot(
  bar_data,
  aes(x = total_pr8_8segment_umi, y = number_of_spots)
) +
  geom_col(width = 0.9, fill = "#2166AC", colour = NA) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_log10(
    breaks = scales::breaks_log(n = 7),
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "Distribution of total PR8 UMI across spatial spots",
    subtitle = paste0(
      "Four samples pooled; exact integer counts for PB2 + PB1 + PA + HA + ",
      "NP + NA + M + NS"
    ),
    x = "Number of all PR8 UMI per spot",
    y = "Number of spots (log10 scale)",
    caption = "Each bar is the number of spots with exactly that total raw PR8 UMI count."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 9, colour = "grey35", hjust = 0)
  )

bar_output_png <- file.path(
  output_dir,
  "PR8_8segment_total_UMI_spot_count_bargraph.png"
)
ggsave(
  bar_output_png, bar_plot,
  width = 9.5, height = 6.5, units = "in", dpi = 300, bg = "white"
)
violin_plot <- ggplot(
  plot_data,
  aes(x = sample_id, y = log1p_total_pr8_8segment_umi, fill = sample_id)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, colour = "grey25") +
  geom_boxplot(
    width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.8,
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "2_2" = "#2166AC", "2_3" = "#67A9CF",
      "4_2" = "#D6604D", "4_3" = "#F4A582"
    ),
    guide = "none"
  ) +
  labs(
    title = "PR8 eight-segment expression per spatial spot",
    subtitle = paste0(
      "PB2 + PB1 + PA + HA + NP + NA + M (M1/M2) + NS (NS1/NEP); ",
      "all spots included"
    ),
    x = "Sample",
    y = "log1p(total PR8 raw UMI per spot)",
    caption = "Blue: Day 4 (2_2, 2_3); orange-red: Day 10 (4_2, 4_3)."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
    axis.text.x = element_text(face = "bold")
  )

output_png <- file.path(output_dir, "PR8_8segment_total_expression_violin.png")
ggsave(
  output_png, violin_plot,
  width = 8.5, height = 6.5, units = "in", dpi = 300, bg = "white"
)

cutoff_data <- data.frame(
  log1p_cutoff = c(0.2, 0.3, 0.5, 0.8, 0.9, 1.0),
  minimum_raw_umi = c(1L, 1L, 1L, 2L, 2L, 2L),
  stringsAsFactors = FALSE
)
cutoff_data$cutoff_label <- sprintf(
  "%.1f  (raw UMI >= %d)",
  cutoff_data$log1p_cutoff,
  cutoff_data$minimum_raw_umi
)
cutoff_colours <- setNames(
  c("#08306B", "#2171B5", "#6BAED6", "#A63603", "#E6550D", "#FDAE6B"),
  cutoff_data$cutoff_label
)

violin_plot_with_cutoffs <- violin_plot +
  geom_hline(
    data = cutoff_data,
    aes(yintercept = log1p_cutoff, colour = cutoff_label),
    linetype = "dashed", linewidth = 0.65, alpha = 0.95,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = cutoff_colours,
    name = "log1p cutoff\n(equivalent raw UMI)"
  ) +
  labs(
    title = "PR8 eight-segment expression with requested cutoffs",
    caption = paste0(
      "Dashed lines show requested log1p cutoffs. Raw UMI is integer-valued: ",
      "cutoffs 0.2/0.3/0.5 all imply raw UMI >= 1; ",
      "0.8/0.9/1.0 all imply raw UMI >= 2."
    )
  ) +
  theme(legend.position = "right")

cutoff_output_png <- file.path(
  output_dir,
  "PR8_8segment_total_expression_violin_with_cutoffs.png"
)
ggsave(
  cutoff_output_png, violin_plot_with_cutoffs,
  width = 10.5, height = 6.5, units = "in", dpi = 300, bg = "white"
)

raw_violin_plot <- ggplot(
  plot_data,
  aes(x = sample_id, y = total_pr8_8segment_umi, fill = sample_id)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, colour = "grey25") +
  geom_boxplot(
    width = 0.13, outlier.shape = NA, fill = "white", alpha = 0.8,
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "2_2" = "#2166AC", "2_3" = "#67A9CF",
      "4_2" = "#D6604D", "4_3" = "#F4A582"
    ),
    guide = "none"
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "PR8 eight-segment raw expression per spatial spot",
    subtitle = paste0(
      "PB2 + PB1 + PA + HA + NP + NA + M (M1/M2) + NS (NS1/NEP); ",
      "all spots included"
    ),
    x = "Sample",
    y = "Total PR8 raw UMI per spot",
    caption = "Blue: Day 4 (2_2, 2_3); orange-red: Day 10 (4_2, 4_3). No log transformation."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
    axis.text.x = element_text(face = "bold")
  )

raw_output_png <- file.path(
  output_dir,
  "PR8_8segment_total_expression_violin_raw_UMI.png"
)
ggsave(
  raw_output_png, raw_violin_plot,
  width = 8.5, height = 6.5, units = "in", dpi = 300, bg = "white"
)

raw_violin_zoom_plot <- raw_violin_plot +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "PR8 eight-segment raw expression per spatial spot (0–100 UMI)",
    caption = paste0(
      "Display zoom: 0–100 raw UMI per spot; values above 100 remain in the ",
      "violin density calculation but are outside the displayed y range."
    )
  )

raw_zoom_output_png <- file.path(
  output_dir,
  "PR8_8segment_total_expression_violin_raw_UMI_zoom_0_100.png"
)
ggsave(
  raw_zoom_output_png, raw_violin_zoom_plot,
  width = 8.5, height = 6.5, units = "in", dpi = 300, bg = "white"
)

raw_violin_zoom_200_plot <- raw_violin_plot +
  coord_cartesian(ylim = c(0, 200)) +
  scale_y_continuous(
    breaks = seq(0, 200, by = 50),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "PR8 eight-segment raw expression per spatial spot (0–200 UMI)",
    caption = paste0(
      "Display zoom: 0–200 raw UMI per spot; values above 200 remain in the ",
      "violin density calculation but are outside the displayed y range."
    )
  )

raw_zoom_200_output_png <- file.path(
  output_dir,
  "PR8_8segment_total_expression_violin_raw_UMI_zoom_0_200.png"
)
ggsave(
  raw_zoom_200_output_png, raw_violin_zoom_200_plot,
  width = 8.5, height = 6.5, units = "in", dpi = 300, bg = "white"
)

writeLines(
  c(
    "PR8 eight-segment spot-expression violin plot",
    "",
    "Eight segments: PB2, PB1, PA, HA, NP, NA, M, and NS.",
    "M is the sum of available M/M1/M2 features; NS is the sum of available NS/NS1/NS2/NEP features.",
    "A segment absent from an input count matrix contributes zero; consult PR8_8segment_feature_mapping_by_sample.tsv.",
    "All spots, including zero-expression spots, are shown using log1p(total raw UMI).",
    "PR8_8segment_total_UMI_spot_count_bargraph.png pools all four samples and counts spots at each exact total raw PR8 UMI value; its y axis uses a log10 display scale."
  ),
  file.path(output_dir, "README_PR8_8segment_violin.txt")
)

message("Saved: ", output_png)
message("Saved: ", cutoff_output_png)
message("Saved: ", raw_output_png)
message("Saved: ", raw_zoom_output_png)
message("Saved: ", raw_zoom_200_output_png)
message("Saved: ", bar_output_png)
