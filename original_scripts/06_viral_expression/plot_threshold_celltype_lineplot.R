#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

get_script_path <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) == 0L) {
    return(file.path(getwd(), "plot_threshold_celltype_lineplot.R"))
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

analysis_dir <- dirname(get_script_path())
input_file <- file.path(
  analysis_dir,
  "thresholded_celltype_barline",
  "Day4_RCTD_celltype_raw_UMI_threshold_series.tsv"
)
output_dir <- file.path(analysis_dir, "threshold_celltype_lineplot")
output_stem <- "Day4_RCTD_infected_percentage_thresholds_01_02_05_10"
plot_width <- 12
plot_height <- 12

selected_thresholds <- c(1L, 2L, 5L, 10L)
threshold_labels <- sprintf("\u2265 %d UMI", selected_thresholds)
names(threshold_labels) <- as.character(selected_thresholds)
line_colours <- setNames(
  c("#BE3035", "#E83C34", "#E3676D", "#E1939C"),
  threshold_labels
)

if (!file.exists(input_file)) {
  stop("Input summary does not exist: ", input_file)
}

summary_table <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "cell_type", "total_spots", "infected_spots", "uninfected_spots",
  "infected_percent", "raw_umi_threshold"
)
missing_columns <- setdiff(required_columns, colnames(summary_table))
if (length(missing_columns) > 0L) {
  stop("Input summary is missing columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(summary_table[c("cell_type", "raw_umi_threshold")])) {
  stop("Input summary has duplicate cell_type/raw_umi_threshold rows.")
}

available_thresholds <- sort(unique(summary_table$raw_umi_threshold))
missing_thresholds <- setdiff(selected_thresholds, available_thresholds)
if (length(missing_thresholds) > 0L) {
  stop("Input summary is missing thresholds: ", paste(missing_thresholds, collapse = ", "))
}

# The threshold-1 block retains the fixed order used by every existing subfolder.
cell_type_order <- summary_table$cell_type[
  summary_table$raw_umi_threshold == selected_thresholds[[1L]]
]
if (length(cell_type_order) == 0L || anyDuplicated(cell_type_order)) {
  stop("Could not determine a unique cell-type order from threshold 1.")
}

selected_data <- summary_table[
  summary_table$raw_umi_threshold %in% selected_thresholds,
  required_columns,
  drop = FALSE
]
selected_data$cell_type <- factor(selected_data$cell_type, levels = cell_type_order)
selected_data$threshold_label <- factor(
  unname(threshold_labels[as.character(selected_data$raw_umi_threshold)]),
  levels = threshold_labels
)
selected_data <- selected_data[
  order(selected_data$raw_umi_threshold, selected_data$cell_type),
  ,
  drop = FALSE
]

expected_rows <- length(cell_type_order) * length(selected_thresholds)
if (nrow(selected_data) != expected_rows || anyNA(selected_data$cell_type)) {
  stop("Selected data do not contain every cell type at all four thresholds.")
}
if (anyNA(selected_data$infected_percent) ||
    any(selected_data$infected_percent < 0 | selected_data$infected_percent > 100)) {
  stop("infected_percent contains missing values or values outside [0, 100].")
}

recalculated_percent <- 100 * selected_data$infected_spots / selected_data$total_spots
if (!isTRUE(all.equal(
  selected_data$infected_percent,
  recalculated_percent,
  tolerance = 1e-10,
  check.attributes = FALSE
))) {
  stop("infected_percent is inconsistent with infected_spots / total_spots.")
}

for (cell_type_name in cell_type_order) {
  cell_values <- selected_data$infected_percent[
    selected_data$cell_type == cell_type_name
  ]
  if (any(diff(cell_values) > 1e-10)) {
    stop("Infected percentage increases with threshold for cell type: ", cell_type_name)
  }
}

base_plot <- ggplot(
  selected_data,
  aes(
    x = cell_type,
    y = infected_percent,
    colour = threshold_label,
    group = threshold_label
  )
) +
  geom_line(linewidth = 1.25, lineend = "round", linejoin = "round") +
  scale_colour_manual(values = line_colours, drop = FALSE) +
  scale_x_discrete(
    drop = FALSE,
    labels = function(labels) {
      vapply(
        labels,
        function(label) paste(strwrap(label, width = 18), collapse = "\n"),
        character(1)
      )
    },
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  labs(
    x = "Cell type",
    y = "Infected spots (%)",
    colour = NULL
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 1.8)
    )
  ) +
  theme_classic(base_size = 13, base_family = "sans") +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    axis.title.x = element_text(size = 14, margin = margin(t = 14)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text.x = element_text(
      angle = 65,
      hjust = 1,
      vjust = 1,
      size = 15,
      lineheight = 0.95,
      colour = "black"
    ),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.55),
    axis.ticks.length = grid::unit(2.5, "mm"),
    axis.line = element_line(colour = "black", linewidth = 0.9),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(18, "mm"),
    legend.key.height = grid::unit(5, "mm"),
    legend.spacing.x = grid::unit(3, "mm"),
    plot.margin = margin(10, 16, 12, 12)
  )

line_plot <- base_plot +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(value) paste0(value, "%"),
    expand = expansion(mult = c(0, 0.015))
  )

if (any(selected_data$infected_percent <= 0)) {
  stop("A log10 y-axis requires all infected_percent values to be greater than zero.")
}
minimum_infected_percent <- min(selected_data$infected_percent)
log10_lower_limit <- min(1, 0.8 * minimum_infected_percent)
log10_plot <- base_plot +
  scale_y_log10(
    limits = c(log10_lower_limit, 100),
    breaks = c(1, 10, 100),
    labels = c("1%", "10%", "100%"),
    expand = expansion(mult = c(0.01, 0.015))
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(output_dir, paste0(output_stem, ".png"))
pdf_file <- file.path(output_dir, paste0(output_stem, ".pdf"))
svg_file <- file.path(output_dir, paste0(output_stem, ".svg"))
tsv_file <- file.path(output_dir, paste0(output_stem, ".tsv"))
log10_output_stem <- paste0(output_stem, "_log10")
log10_png_file <- file.path(output_dir, paste0(log10_output_stem, ".png"))
log10_pdf_file <- file.path(output_dir, paste0(log10_output_stem, ".pdf"))
log10_svg_file <- file.path(output_dir, paste0(log10_output_stem, ".svg"))

ggsave(
  png_file,
  plot = line_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  pdf_file,
  plot = line_plot,
  device = grDevices::cairo_pdf,
  width = plot_width,
  height = plot_height,
  units = "in",
  bg = "white"
)
ggsave(
  svg_file,
  plot = line_plot,
  device = grDevices::svg,
  width = plot_width,
  height = plot_height,
  units = "in",
  bg = "white"
)
ggsave(
  log10_png_file,
  plot = log10_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  log10_pdf_file,
  plot = log10_plot,
  device = grDevices::cairo_pdf,
  width = plot_width,
  height = plot_height,
  units = "in",
  bg = "white"
)
ggsave(
  log10_svg_file,
  plot = log10_plot,
  device = grDevices::svg,
  width = plot_width,
  height = plot_height,
  units = "in",
  bg = "white"
)

export_data <- selected_data
export_data$cell_type <- as.character(export_data$cell_type)
export_data$threshold_label <- as.character(export_data$threshold_label)
write.table(
  export_data,
  tsv_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  c(
    "Day 4 RCTD cell-type infected-percentage line plot",
    "",
    paste0("Input: ", input_file),
    "Samples: pooled Day 4 samples 2_2 and 2_3.",
    paste0(
      "Included RCTD spot classes: singlet, doublet_uncertain, ",
      "doublet_certain; cell type is first_type."
    ),
    paste0(
      "Definition: infected = total PR8 raw UMI >= threshold; ",
      "infected percentage = 100 * infected spots / total spots."
    ),
    "Displayed thresholds: 1, 2, 5, and 10 raw UMI.",
    paste0("Figure dimensions: ", plot_width, " x ", plot_height, " inches (1:1)."),
    paste0("Cell-type order: ", paste(cell_type_order, collapse = " | ")),
    "The figure has no title, subtitle, caption, point labels, or other annotations.",
    paste0("PNG: ", basename(png_file)),
    paste0("PDF (vector): ", basename(pdf_file)),
    paste0("SVG (vector): ", basename(svg_file)),
    paste0(
      "Log10 y-axis: ticks at 1%, 10%, and 100%; lower plotting limit = ",
      format(log10_lower_limit, digits = 4),
      "% to retain the observed minimum of ",
      format(minimum_infected_percent, digits = 4),
      "%."
    ),
    paste0("Log10 PNG: ", basename(log10_png_file)),
    paste0("Log10 PDF (vector): ", basename(log10_pdf_file)),
    paste0("Log10 SVG (vector): ", basename(log10_svg_file)),
    paste0("Plotted data: ", basename(tsv_file)),
    paste0("Generated by: ", basename(get_script_path()))
  ),
  file.path(output_dir, "README.txt")
)

message("Saved: ", png_file)
message("Saved: ", pdf_file)
message("Saved: ", svg_file)
message("Saved: ", log10_png_file)
message("Saved: ", log10_pdf_file)
message("Saved: ", log10_svg_file)
message("Saved: ", tsv_file)
