#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_dir <- "/data/projects/punim2394/lung_curio_mouse_pr8"
output_dir <- file.path(project_dir, "viralgeneassess")
input_file <- file.path(output_dir, "Day4_RCTD_celltype_infection_pooled.tsv")
output_file <- file.path(
  output_dir,
  "Day4_RCTD_celltype_infected_bar_with_percentage_line.png"
)

if (!file.exists(input_file)) {
  stop("Input summary is missing: ", input_file)
}

plot_data <- read.delim(input_file, check.names = FALSE)
required_columns <- c("cell_type", "infected_spots", "infected_percent")
missing_columns <- setdiff(required_columns, colnames(plot_data))
if (length(missing_columns) > 0L) {
  stop("Input summary is missing columns: ", paste(missing_columns, collapse = ", "))
}

plot_data$cell_type <- factor(plot_data$cell_type, levels = plot_data$cell_type)

percentage_plot <- ggplot(
  plot_data,
  aes(x = cell_type, y = infected_percent, group = 1)
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
    limits = c(55, 92),
    breaks = seq(60, 90, by = 10),
    labels = function(value) paste0(value, "%"),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(x = NULL, y = "Infected spots (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(8, 15, 0, 12)
  )

bar_plot <- ggplot(
  plot_data,
  aes(x = cell_type, y = infected_spots)
) +
  geom_col(
    width = 0.68,
    fill = "#D73027", colour = "grey25", linewidth = 0.3
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(x = NULL, y = "Number of infected spots") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(0, 15, 8, 12)
  )

combined_plot <- percentage_plot / bar_plot +
  plot_layout(heights = c(1.25, 3.0)) +
  plot_annotation(
    title = "Day 4 PR8 infection by RCTD cell type",
    subtitle = paste0(
      "Replicates 2_2 and 2_3 pooled | accepted RCTD spots only | ",
      "infected = total PR8 UMI > 0"
    ),
    caption = paste0(
      "Bars show infected spot counts. The aligned line and labels show the ",
      "infected percentage within each cell type. Cell type is the dominant ",
      "RCTD label (first_type)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9, colour = "grey35", hjust = 0),
      plot.margin = margin(12, 12, 12, 12)
    )
  )

ggsave(
  output_file,
  combined_plot,
  width = 13,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Saved: ", output_file)
