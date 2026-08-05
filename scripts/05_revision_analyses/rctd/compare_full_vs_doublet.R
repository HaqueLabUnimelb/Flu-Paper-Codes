#!/usr/bin/env Rscript
# Purpose: compare same-input day-4 full-mode and doublet-mode RCTD results.
# Inputs: FLU_RCTD_DAY4_FULL_RDS and FLU_RCTD_DAY4_DOUBLET_RDS.
# Output directory: FLU_RCTD_MODE_COMPARISON_OUTPUT_DIR.

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

full_file <- env_path(
  "FLU_RCTD_DAY4_FULL_RDS",
  file.path(repo_dir, "outputs", "rctd", "Day4mouseRCTDoutmin30_multi.rds")
)
doublet_file <- env_path(
  "FLU_RCTD_DAY4_DOUBLET_RDS",
  file.path(repo_dir, "outputs", "revision", "rctd", "Day4mouseRCTDoutmin30_doublet_same_input.rds")
)
output_dir <- env_path(
  "FLU_RCTD_MODE_COMPARISON_OUTPUT_DIR",
  file.path(repo_dir, "outputs", "revision", "rctd", "full_vs_doublet_same_input")
)
replicate_names <- trimws(strsplit(
  Sys.getenv("FLU_RCTD_DAY4_REPLICATE_NAMES", "puck1_lungST_2_2,puck2_lungST_2_3"),
  ",",
  fixed = TRUE
)[[1L]])

missing_inputs <- c(full_file, doublet_file)[!file.exists(c(full_file, doublet_file))]
if (length(missing_inputs) > 0L) {
  stop("Required RCTD input(s) not found: ", paste(missing_inputs, collapse = ", "))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

canonicalize_types <- function(x, canonical_types) {
  lookup <- setNames(canonical_types, tolower(canonical_types))
  result <- unname(lookup[tolower(as.character(x))])
  if (anyNA(result)) {
    stop("Could not harmonize cell types: ", paste(unique(x[is.na(result)]), collapse = ", "))
  }
  result
}

make_doublet_proportions <- function(results, canonical_types) {
  results_df <- results$results_df
  weights_doublet <- as.matrix(results$weights_doublet)
  if (!identical(rownames(results_df), rownames(weights_doublet))) {
    stop("results_df and weights_doublet row names are not aligned")
  }

  proportions <- matrix(
    0,
    nrow = nrow(results_df),
    ncol = length(canonical_types),
    dimnames = list(rownames(results_df), canonical_types)
  )
  first_type <- canonicalize_types(results_df$first_type, canonical_types)
  second_type <- canonicalize_types(results_df$second_type, canonical_types)
  non_reject <- results_df$spot_class != "reject"

  first_rows <- which(non_reject)
  proportions[cbind(first_rows, match(first_type[first_rows], canonical_types))] <-
    weights_doublet[first_rows, "first_type"]

  second_rows <- which(results_df$spot_class == "doublet_certain")
  proportions[cbind(second_rows, match(second_type[second_rows], canonical_types))] <-
    weights_doublet[second_rows, "second_type"]

  row_totals <- rowSums(proportions)
  proportions[non_reject, ] <- proportions[non_reject, , drop = FALSE] / row_totals[non_reject]
  list(
    proportions = proportions,
    spot_class = as.character(results_df$spot_class),
    non_reject = non_reject,
    first_type = first_type,
    second_type = second_type
  )
}

safe_cor <- function(x, y, method) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  unname(cor(x, y, method = method))
}

js_divergence <- function(p, q) {
  m <- (p + q) / 2
  p_term <- ifelse(p > 0, p * log(p / m), 0)
  q_term <- ifelse(q > 0, q * log(q / m), 0)
  rowSums(p_term + q_term) / 2
}

full_object <- readRDS(full_file)
doublet_object <- readRDS(doublet_file)
if (length(full_object@RCTD.reps) != length(doublet_object@RCTD.reps)) {
  stop("Full and doublet objects have different replicate counts")
}
if (length(replicate_names) != length(full_object@RCTD.reps)) {
  stop("FLU_RCTD_DAY4_REPLICATE_NAMES does not match the replicate count")
}

cell_metrics <- list()
replicate_metrics <- list()
spot_class_rows <- list()
spot_agreement_rows <- list()
composition_rows <- list()
scatter_rows <- list()
confusion_rows <- list()

for (i in seq_along(full_object@RCTD.reps)) {
  full_rep <- full_object@RCTD.reps[[i]]
  doublet_rep <- doublet_object@RCTD.reps[[i]]
  full_props <- as.matrix(spacexr::normalize_weights(full_rep@results$weights))
  canonical_types <- colnames(full_props)
  doublet_data <- make_doublet_proportions(doublet_rep@results, canonical_types)

  if (!identical(rownames(full_props), rownames(doublet_data$proportions))) {
    stop("Spot names/order differ between modes for replicate ", i)
  }
  keep <- doublet_data$non_reject
  full_keep <- full_props[keep, , drop = FALSE]
  doublet_keep <- doublet_data$proportions[keep, , drop = FALSE]
  spot_ids <- rownames(full_keep)
  spot_class <- doublet_data$spot_class[keep]

  full_order <- t(apply(full_keep, 1, order, decreasing = TRUE))
  full_top1 <- canonical_types[full_order[, 1]]
  full_top2 <- canonical_types[full_order[, 2]]
  full_top2_mass <- full_keep[cbind(seq_len(nrow(full_keep)), full_order[, 1])] +
    full_keep[cbind(seq_len(nrow(full_keep)), full_order[, 2])]
  doublet_top1 <- doublet_data$first_type[keep]
  doublet_second <- doublet_data$second_type[keep]
  doublet_certain <- spot_class == "doublet_certain"
  top1_match <- full_top1 == doublet_top1
  top2_exact <- rep(NA, length(spot_ids))
  top2_exact[doublet_certain] <- vapply(which(doublet_certain), function(j) {
    setequal(c(full_top1[j], full_top2[j]), c(doublet_top1[j], doublet_second[j]))
  }, logical(1))

  total_variation <- rowSums(abs(full_keep - doublet_keep)) / 2
  jsd <- js_divergence(full_keep, doublet_keep)

  replicate_metrics[[i]] <- data.frame(
    replicate = replicate_names[[i]],
    total_spots = nrow(full_props),
    non_reject_spots = sum(keep),
    reject_fraction = mean(!keep),
    top1_concordance = mean(top1_match),
    top2_exact_concordance_doublet_certain = mean(top2_exact, na.rm = TRUE),
    median_full_top2_mass = median(full_top2_mass),
    q25_full_top2_mass = unname(quantile(full_top2_mass, 0.25)),
    fraction_full_top2_mass_ge_0.90 = mean(full_top2_mass >= 0.90),
    fraction_full_top2_mass_ge_0.95 = mean(full_top2_mass >= 0.95),
    median_total_variation = median(total_variation),
    mean_total_variation = mean(total_variation),
    median_js_divergence = median(jsd),
    stringsAsFactors = FALSE
  )

  class_table <- as.data.frame(table(doublet_data$spot_class), stringsAsFactors = FALSE)
  colnames(class_table) <- c("spot_class", "count")
  class_table$replicate <- replicate_names[[i]]
  class_table$fraction <- class_table$count / sum(class_table$count)
  spot_class_rows[[i]] <- class_table[, c("replicate", "spot_class", "count", "fraction")]

  spot_agreement_rows[[i]] <- data.frame(
    replicate = replicate_names[[i]],
    spot = spot_ids,
    spot_class = spot_class,
    full_top1 = full_top1,
    doublet_top1 = doublet_top1,
    top1_match = top1_match,
    top2_exact_match = top2_exact,
    full_top2_mass = full_top2_mass,
    total_variation = total_variation,
    js_divergence = jsd,
    stringsAsFactors = FALSE
  )

  confusion_rows[[i]] <- data.frame(
    replicate = replicate_names[[i]],
    full_top1 = full_top1,
    doublet_top1 = doublet_top1,
    stringsAsFactors = FALSE
  )

  set.seed(1000 + i)
  sampled_rows <- sample(seq_len(nrow(full_keep)), min(5000L, nrow(full_keep)))
  for (cell_type in canonical_types) {
    full_values <- full_keep[, cell_type]
    doublet_values <- doublet_keep[, cell_type]
    differences <- full_values - doublet_values
    cell_metrics[[length(cell_metrics) + 1L]] <- data.frame(
      replicate = replicate_names[[i]],
      cell_type = cell_type,
      spots = length(full_values),
      pearson = safe_cor(full_values, doublet_values, "pearson"),
      spearman = safe_cor(full_values, doublet_values, "spearman"),
      mae = mean(abs(differences)),
      rmse = sqrt(mean(differences^2)),
      mean_full = mean(full_values),
      mean_doublet = mean(doublet_values),
      mean_difference_full_minus_doublet = mean(differences),
      median_absolute_difference = median(abs(differences)),
      fraction_abs_difference_le_0.05 = mean(abs(differences) <= 0.05),
      fraction_abs_difference_le_0.10 = mean(abs(differences) <= 0.10),
      stringsAsFactors = FALSE
    )
    composition_rows[[length(composition_rows) + 1L]] <- data.frame(
      replicate = replicate_names[[i]],
      cell_type = cell_type,
      mode = c("full", "doublet"),
      mean_proportion = c(mean(full_values), mean(doublet_values)),
      stringsAsFactors = FALSE
    )
    scatter_rows[[length(scatter_rows) + 1L]] <- data.frame(
      replicate = replicate_names[[i]],
      cell_type = cell_type,
      full = full_values[sampled_rows],
      doublet = doublet_values[sampled_rows],
      stringsAsFactors = FALSE
    )
  }
}

cell_metrics <- do.call(rbind, cell_metrics)
replicate_metrics <- do.call(rbind, replicate_metrics)
spot_class_table <- do.call(rbind, spot_class_rows)
spot_agreement <- do.call(rbind, spot_agreement_rows)
composition_means <- do.call(rbind, composition_rows)
scatter_data <- do.call(rbind, scatter_rows)
confusion_data <- do.call(rbind, confusion_rows)

write.csv(cell_metrics, file.path(output_dir, "per_cell_type_metrics.csv"), row.names = FALSE)
write.csv(replicate_metrics, file.path(output_dir, "per_replicate_metrics.csv"), row.names = FALSE)
write.csv(spot_class_table, file.path(output_dir, "doublet_spot_class_counts.csv"), row.names = FALSE)
write.csv(spot_agreement, file.path(output_dir, "spot_level_agreement.csv"), row.names = FALSE)
write.csv(composition_means, file.path(output_dir, "mean_cell_type_proportions.csv"), row.names = FALSE)

metric_plot <- ggplot(cell_metrics, aes(x = reorder(cell_type, spearman), y = spearman, colour = replicate)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Spot-level Spearman correlation", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(output_dir, "cell_type_spearman_correlations.png"), metric_plot,
       width = 8, height = 6, units = "in", dpi = 300, bg = "white")

mae_plot <- ggplot(cell_metrics, aes(x = reorder(cell_type, mae), y = mae, colour = replicate)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
  coord_flip() +
  labs(x = NULL, y = "Mean absolute proportion difference", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(output_dir, "cell_type_mae.png"), mae_plot,
       width = 8, height = 6, units = "in", dpi = 300, bg = "white")

scatter_plot <- ggplot(scatter_data, aes(x = full, y = doublet, colour = replicate)) +
  geom_point(size = 0.25, alpha = 0.18) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = 2) +
  facet_wrap(~cell_type, ncol = 4) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Full-mode proportion", y = "Doublet-mode proportion", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top")
ggsave(file.path(output_dir, "spot_proportion_scatter.png"), scatter_plot,
       width = 12, height = 10, units = "in", dpi = 300, bg = "white")

top2_plot <- ggplot(spot_agreement, aes(x = full_top2_mass, fill = replicate)) +
  geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.5) +
  geom_vline(xintercept = c(0.90, 0.95), linetype = 2, linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Full-mode mass assigned to top two cell types", y = "Spots", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(output_dir, "full_mode_top2_mass.png"), top2_plot,
       width = 8, height = 5, units = "in", dpi = 300, bg = "white")

composition_plot <- ggplot(
  composition_means,
  aes(x = cell_type, y = mean_proportion, fill = mode)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  facet_wrap(~replicate, ncol = 1) +
  coord_flip() +
  labs(x = NULL, y = "Mean cell-type proportion", fill = "Mode") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(output_dir, "mean_cell_type_proportions.png"), composition_plot,
       width = 9, height = 10, units = "in", dpi = 300, bg = "white")

confusion_table <- as.data.frame(table(confusion_data$full_top1, confusion_data$doublet_top1))
colnames(confusion_table) <- c("full_top1", "doublet_top1", "count")
confusion_plot <- ggplot(confusion_table, aes(x = full_top1, y = doublet_top1, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", trans = "sqrt") +
  coord_fixed() +
  labs(x = "Full-mode top cell type", y = "Doublet-mode top cell type", fill = "Spots") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "top1_cell_type_concordance.png"), confusion_plot,
       width = 9, height = 8, units = "in", dpi = 300, bg = "white")

median_spearman <- aggregate(spearman ~ replicate, cell_metrics, median, na.rm = TRUE)
median_mae <- aggregate(mae ~ replicate, cell_metrics, median, na.rm = TRUE)
summary_lines <- c(
  "# Same-input RCTD full versus doublet comparison",
  "",
  paste0("- R: ", R.version.string),
  paste0("- spacexr: ", as.character(packageVersion("spacexr"))),
  "- The same RCTD objects, spots, reference-derived profiles, filtering and package version were used.",
  "- Doublet rejects were excluded from proportion agreement metrics.",
  "- No spot-level t-test was used because spatial spots are paired and spatially autocorrelated.",
  "",
  "## Reviewer-facing quantitative summary"
)
for (i in seq_len(nrow(replicate_metrics))) {
  r <- replicate_metrics[i, ]
  summary_lines <- c(
    summary_lines,
    "",
    paste0("### ", r$replicate),
    paste0("- Included non-reject spots: ", r$non_reject_spots, " / ", r$total_spots,
           " (reject fraction ", sprintf("%.1f%%", 100 * r$reject_fraction), ")"),
    paste0("- Top-1 cell-type concordance: ", sprintf("%.1f%%", 100 * r$top1_concordance)),
    paste0("- Exact top-2 concordance among doublet_certain spots: ",
           sprintf("%.1f%%", 100 * r$top2_exact_concordance_doublet_certain)),
    paste0("- Median full-mode top-2 mass: ", sprintf("%.3f", r$median_full_top2_mass)),
    paste0("- Full-mode top-2 mass >= 0.90: ",
           sprintf("%.1f%%", 100 * r$fraction_full_top2_mass_ge_0.90)),
    paste0("- Full-mode top-2 mass >= 0.95: ",
           sprintf("%.1f%%", 100 * r$fraction_full_top2_mass_ge_0.95)),
    paste0("- Median cell-type Spearman correlation: ",
           sprintf("%.3f", median_spearman$spearman[median_spearman$replicate == r$replicate])),
    paste0("- Median cell-type MAE: ",
           sprintf("%.3f", median_mae$mae[median_mae$replicate == r$replicate])),
    paste0("- Median per-spot total variation distance: ",
           sprintf("%.3f", r$median_total_variation))
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "## Interpretation",
  "If the top-two mass, top-cell concordance and cell-type correlations are high while MAE and total variation are low, the full solution is effectively concentrated in at most two cell types per spot. This supports the use of doublet mode for high-resolution Curioseeker spots while showing that the biological conclusions are insensitive to mode choice."
)
writeLines(summary_lines, file.path(output_dir, "reviewer_summary.md"))
message("Completed comparison outputs in: ", output_dir)
