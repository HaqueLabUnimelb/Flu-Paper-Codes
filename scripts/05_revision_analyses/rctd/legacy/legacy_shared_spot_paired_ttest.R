#!/usr/bin/env Rscript
# Legacy exploratory comparison. The source objects were not generated from
# identical inputs, and spot-level tests do not provide biological-replicate
# inference. Prefer the same-input comparison in the parent directory.

suppressPackageStartupMessages(library(spacexr))

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Unable to determine the script location. Run this file with Rscript.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
repo_dir <- normalizePath(file.path(script_dir, "..", "..", "..", ".."), mustWork = TRUE)

env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) {
    stop("Required path is empty: ", name)
  }
  value
}

doublet_file <- env_path(
  "FLU_RCTD_LEGACY_DOUBLET_RDS",
  file.path(repo_dir, "data", "processed", "Day4mouseRCTDoutmin30_editv1")
)
full_file <- env_path(
  "FLU_RCTD_DAY4_FULL_RDS",
  file.path(repo_dir, "outputs", "rctd", "Day4mouseRCTDoutmin30_multi.rds")
)
output_dir <- env_path(
  "FLU_RCTD_LEGACY_COMPARISON_OUTPUT_DIR",
  file.path(repo_dir, "outputs", "revision", "rctd", "legacy_shared_spot_paired_ttest")
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
    stop("results_df and weights_doublet rows are not aligned")
  }

  proportions <- matrix(
    0,
    nrow = nrow(results_df),
    ncol = length(canonical_types),
    dimnames = list(rownames(results_df), canonical_types)
  )
  first_type <- canonicalize_types(results_df$first_type, canonical_types)
  second_type <- canonicalize_types(results_df$second_type, canonical_types)
  non_reject <- as.character(results_df$spot_class) != "reject"

  first_rows <- which(non_reject)
  proportions[cbind(first_rows, match(first_type[first_rows], canonical_types))] <-
    weights_doublet[first_rows, "first_type"]

  second_rows <- which(as.character(results_df$spot_class) == "doublet_certain")
  proportions[cbind(second_rows, match(second_type[second_rows], canonical_types))] <-
    weights_doublet[second_rows, "second_type"]

  row_totals <- rowSums(proportions)
  proportions[non_reject, ] <- proportions[non_reject, , drop = FALSE] / row_totals[non_reject]

  list(
    proportions = proportions,
    spot_class = as.character(results_df$spot_class),
    non_reject = non_reject
  )
}

run_paired_tests <- function(full_table, doublet_table, cell_types, scope) {
  rows <- lapply(cell_types, function(cell_type) {
    full_values <- full_table[[cell_type]]
    doublet_values <- doublet_table[[cell_type]]
    differences <- full_values - doublet_values
    test <- t.test(full_values, doublet_values, paired = TRUE)
    data.frame(
      scope = scope,
      cell_type = cell_type,
      paired_spots = length(differences),
      mean_full = mean(full_values),
      mean_doublet = mean(doublet_values),
      mean_difference_full_minus_doublet = mean(differences),
      median_difference_full_minus_doublet = median(differences),
      mean_absolute_difference = mean(abs(differences)),
      t_statistic = unname(test$statistic),
      degrees_of_freedom = unname(test$parameter),
      confidence_interval_low = test$conf.int[[1]],
      confidence_interval_high = test$conf.int[[2]],
      p_value = test$p.value,
      cohen_dz = mean(differences) / sd(differences),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$p_adjust_bh <- p.adjust(result$p_value, method = "BH")
  result
}

doublet_object <- readRDS(doublet_file)
full_object <- readRDS(full_file)
if (length(doublet_object@RCTD.reps) != length(full_object@RCTD.reps)) {
  stop("The objects contain different numbers of replicates")
}
if (length(replicate_names) != length(full_object@RCTD.reps)) {
  stop("FLU_RCTD_DAY4_REPLICATE_NAMES does not match the replicate count")
}

full_tables <- list()
doublet_tables <- list()
test_results <- list()
overlap_rows <- list()

for (i in seq_along(full_object@RCTD.reps)) {
  full_weights <- as.matrix(
    spacexr::normalize_weights(full_object@RCTD.reps[[i]]@results$weights)
  )
  canonical_types <- colnames(full_weights)
  doublet_data <- make_doublet_proportions(
    doublet_object@RCTD.reps[[i]]@results,
    canonical_types
  )

  shared_spots <- intersect(rownames(doublet_data$proportions), rownames(full_weights))
  shared_non_reject <- shared_spots[
    doublet_data$spot_class[match(shared_spots, rownames(doublet_data$proportions))] != "reject"
  ]
  if (length(shared_non_reject) == 0L) {
    stop("No shared non-reject spots for replicate ", i)
  }

  full_props <- full_weights[shared_non_reject, canonical_types, drop = FALSE]
  doublet_props <- doublet_data$proportions[shared_non_reject, canonical_types, drop = FALSE]
  spot_class <- doublet_data$spot_class[
    match(shared_non_reject, rownames(doublet_data$proportions))
  ]

  full_table <- data.frame(
    replicate = replicate_names[[i]],
    spot = shared_non_reject,
    doublet_spot_class = spot_class,
    full_props,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  doublet_table <- data.frame(
    replicate = replicate_names[[i]],
    spot = shared_non_reject,
    doublet_spot_class = spot_class,
    doublet_props,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!identical(full_table[, c("replicate", "spot")], doublet_table[, c("replicate", "spot")])) {
    stop("Paired table row identifiers are not identical for replicate ", i)
  }
  full_tables[[i]] <- full_table
  doublet_tables[[i]] <- doublet_table
  test_results[[i]] <- run_paired_tests(
    full_table,
    doublet_table,
    canonical_types,
    replicate_names[[i]]
  )
  overlap_rows[[i]] <- data.frame(
    replicate = replicate_names[[i]],
    old_doublet_spots = nrow(doublet_data$proportions),
    new_full_spots = nrow(full_weights),
    shared_spots = length(shared_spots),
    shared_non_reject_spots_used = length(shared_non_reject),
    stringsAsFactors = FALSE
  )
}

full_table <- do.call(rbind, full_tables)
doublet_table <- do.call(rbind, doublet_tables)
cell_types <- colnames(full_table)[-(1:3)]

if (!identical(full_table[, c("replicate", "spot")], doublet_table[, c("replicate", "spot")])) {
  stop("Final paired table row identifiers are not identical")
}
if (!all(abs(rowSums(full_table[, cell_types]) - 1) < 1e-8)) {
  stop("Full-mode proportions do not sum to one")
}
if (!all(abs(rowSums(doublet_table[, cell_types]) - 1) < 1e-8)) {
  stop("Doublet-mode proportions do not sum to one")
}

test_results[[length(test_results) + 1L]] <- run_paired_tests(
  full_table,
  doublet_table,
  cell_types,
  "combined_both_pucks"
)
t_test_results <- do.call(rbind, test_results)
overlap_summary <- do.call(rbind, overlap_rows)

write.csv(
  full_table,
  file.path(output_dir, "full_mode_shared_spot_proportions.csv"),
  row.names = FALSE
)
write.csv(
  doublet_table,
  file.path(output_dir, "doublet_mode_shared_spot_proportions.csv"),
  row.names = FALSE
)
write.csv(
  t_test_results,
  file.path(output_dir, "paired_t_test_results.csv"),
  row.names = FALSE
)
write.csv(
  overlap_summary,
  file.path(output_dir, "spot_overlap_summary.csv"),
  row.names = FALSE
)

note <- c(
  "Exploratory legacy shared-spot comparison.",
  "The two source objects were not produced from identical spot sets/preprocessing, so differences cannot be attributed solely to RCTD mode.",
  "Only shared spots with a non-reject doublet classification were included.",
  "The two proportion CSV files have identical row order and are the exact inputs to the paired t-tests.",
  "Spot-level p-values should not be treated as biological-replicate inference because spatial spots are autocorrelated and only two pucks are available.",
  "Use effect sizes and the separately submitted same-input comparison as the primary reviewer evidence."
)
writeLines(note, file.path(output_dir, "IMPORTANT_README.txt"))

message("Full table rows: ", nrow(full_table))
message("Doublet table rows: ", nrow(doublet_table))
message("Completed exploratory paired t-tests in: ", output_dir)
