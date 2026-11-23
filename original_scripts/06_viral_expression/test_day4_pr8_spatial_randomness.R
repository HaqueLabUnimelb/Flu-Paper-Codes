#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(RANN)
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
    c("annotated_lungST_2_2_QCed.rds", "annotated_lungST_2_3_QCed.rds")
  ),
  stringsAsFactors = FALSE
)

get_option <- function(name, default, cast) {
  prefix <- paste0("--", name, "=")
  hits <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hits) == 0L) {
    return(default)
  }
  if (length(hits) > 1L) {
    stop("Option supplied more than once: --", name)
  }
  cast(sub(prefix, "", hits, fixed = TRUE))
}

n_permutations <- get_option("n-permutations", 1999L, as.integer)
k_neighbors <- get_option("k-neighbors", 6L, as.integer)
library_strata_count <- get_option("library-strata", 10L, as.integer)
minimum_group_size <- get_option("minimum-group-size", 100L, as.integer)
random_seed <- get_option("seed", 20260803L, as.integer)

if (is.na(n_permutations) || n_permutations < 99L) {
  stop("--n-permutations must be at least 99.")
}
if (is.na(k_neighbors) || k_neighbors < 1L) {
  stop("--k-neighbors must be a positive integer.")
}
if (is.na(library_strata_count) || library_strata_count < 1L) {
  stop("--library-strata must be a positive integer.")
}
if (is.na(minimum_group_size) || minimum_group_size < 2L) {
  stop("--minimum-group-size must be at least 2.")
}
stopifnot(all(file.exists(samples$input_file)))

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
  counts <- LayerData(object, assay = assay_name, layer = count_layers)
  list(matrix = counts, assay = assay_name, layer = count_layers)
}

detect_pr8_features <- function(feature_names) {
  tagged <- grepl(
    "(^|[-_.:|/])pr8([-_.:|/]|$)|puerto[ _/-]*rico|influenza[ _-]*a|h1n1",
    feature_names, ignore.case = TRUE, perl = TRUE
  ) | grepl(
    "(^|[-_.:|])iav([-_.:|]|$)",
    feature_names, ignore.case = TRUE, perl = TRUE
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
    stop("No spatial reduction found in object.")
  }
  coordinates <- Embeddings(object, reduction = spatial_reduction[[1L]])
  if (ncol(coordinates) < 2L) {
    stop("Spatial reduction has fewer than two dimensions.")
  }
  coordinates[, 1:2, drop = FALSE]
}

make_knn_edges <- function(coordinates, k) {
  n <- nrow(coordinates)
  if (n <= k) {
    stop("Number of spots must exceed k-neighbors.")
  }
  nearest <- RANN::nn2(coordinates, k = k + 1L)$nn.idx[, -1L, drop = FALSE]
  from <- rep(seq_len(n), each = k)
  to <- as.vector(t(nearest))
  edge_a <- pmin(from, to)
  edge_b <- pmax(from, to)
  non_self <- edge_a != edge_b
  edge_a <- edge_a[non_self]
  edge_b <- edge_b[non_self]
  key <- unique(edge_a + (edge_b - 1) * n)
  data.frame(
    from = as.integer(((key - 1) %% n) + 1),
    to = as.integer(floor((key - 1) / n) + 1),
    stringsAsFactors = FALSE
  )
}

moran_i <- function(values, edges) {
  centered <- values - mean(values)
  denominator <- sum(centered^2)
  if (denominator == 0) {
    return(NA_real_)
  }
  length(values) * sum(centered[edges$from] * centered[edges$to]) /
    (nrow(edges) * denominator)
}

bb_join_counts <- function(values, edges, thresholds) {
  edge_minimum <- pmin(values[edges$from], values[edges$to])
  maximum <- max(values)
  frequency <- tabulate(as.integer(edge_minimum) + 1L, nbins = maximum + 1L)
  at_least <- rev(cumsum(rev(frequency)))
  as.numeric(at_least[as.integer(thresholds) + 1L])
}

make_rank_strata <- function(values, number_of_strata) {
  n <- length(values)
  strata <- ceiling(rank(values, ties.method = "average") * number_of_strata / n)
  pmax(1L, pmin(number_of_strata, as.integer(strata)))
}

permute_within_groups <- function(values, index_groups) {
  result <- values
  for (indices in index_groups) {
    result[indices] <- sample(values[indices], length(indices), replace = FALSE)
  }
  result
}

analyse_sample <- function(sample_id, input_file, sample_index) {
  message("Reading sample ", sample_id, ": ", input_file)
  object <- readRDS(input_file)
  counts_info <- get_counts(object)
  counts <- counts_info$matrix
  pr8_features <- detect_pr8_features(rownames(counts))
  if (length(pr8_features) == 0L) {
    stop("No PR8 features detected in sample ", sample_id, ".")
  }

  viral_umi <- as.numeric(Matrix::colSums(counts[pr8_features, , drop = FALSE]))
  names(viral_umi) <- colnames(counts)
  total_rna_umi <- as.numeric(Matrix::colSums(counts))
  names(total_rna_umi) <- colnames(counts)
  coordinates <- get_spatial_coordinates(object)
  if (!all(names(viral_umi) %in% rownames(coordinates))) {
    stop("Some RNA spots have no spatial coordinates in sample ", sample_id, ".")
  }
  coordinates <- coordinates[names(viral_umi), , drop = FALSE]

  edges <- make_knn_edges(coordinates, k_neighbors)
  positive_values <- sort(unique(viral_umi[viral_umi > 0]))
  positive_count <- vapply(
    positive_values, function(threshold) sum(viral_umi >= threshold), integer(1L)
  )
  negative_count <- length(viral_umi) - positive_count
  keep <- positive_count >= minimum_group_size & negative_count >= minimum_group_size
  thresholds <- positive_values[keep]
  positive_count <- positive_count[keep]
  negative_count <- negative_count[keep]
  if (length(thresholds) == 0L) {
    stop("No eligible thresholds remain in sample ", sample_id, ".")
  }

  observed_joins <- bb_join_counts(viral_umi, edges, thresholds)
  observed_moran <- moran_i(log1p(viral_umi), edges)

  strata <- make_rank_strata(total_rna_umi, library_strata_count)
  index_groups <- split(seq_along(strata), strata)
  permutation_joins <- matrix(
    NA_real_, nrow = n_permutations, ncol = length(thresholds)
  )
  permutation_moran <- numeric(n_permutations)

  set.seed(random_seed + sample_index - 1L)
  for (permutation_index in seq_len(n_permutations)) {
    permuted_umi <- permute_within_groups(viral_umi, index_groups)
    permutation_joins[permutation_index, ] <-
      bb_join_counts(permuted_umi, edges, thresholds)
    permutation_moran[[permutation_index]] <-
      moran_i(log1p(permuted_umi), edges)
    if (permutation_index %% 100L == 0L) {
      message(
        sample_id, ": completed ", permutation_index, "/",
        n_permutations, " permutations"
      )
    }
  }

  null_mean <- colMeans(permutation_joins)
  null_sd <- apply(permutation_joins, 2L, sd)
  observed_z <- (observed_joins - null_mean) / null_sd
  permutation_z <- sweep(permutation_joins, 2L, null_mean, "-")
  permutation_z <- sweep(permutation_z, 2L, null_sd, "/")
  permutation_z[, !is.finite(observed_z)] <- NA_real_
  max_null_z <- apply(permutation_z, 1L, max, na.rm = TRUE)

  raw_p <- vapply(
    seq_along(thresholds),
    function(index) {
      (1 + sum(permutation_joins[, index] >= observed_joins[[index]])) /
        (n_permutations + 1)
    },
    numeric(1L)
  )
  fwer_p <- vapply(
    observed_z,
    function(z_value) {
      if (!is.finite(z_value)) {
        return(NA_real_)
      }
      (1 + sum(max_null_z >= z_value)) / (n_permutations + 1)
    },
    numeric(1L)
  )

  threshold_results <- data.frame(
    sample_id = sample_id,
    threshold_raw_umi = thresholds,
    threshold_log1p = log1p(thresholds),
    positive_spots = positive_count,
    negative_spots = negative_count,
    positive_percent = 100 * positive_count / length(viral_umi),
    observed_positive_positive_edges = observed_joins,
    null_mean_positive_positive_edges = null_mean,
    observed_to_null_ratio = observed_joins / null_mean,
    join_count_z = observed_z,
    permutation_p_unadjusted = raw_p,
    permutation_p_maxT_FWER = fwer_p,
    significant_clustered_FWER_0.05 =
      is.finite(observed_z) & observed_z > 0 & fwer_p <= 0.05,
    stringsAsFactors = FALSE
  )

  moran_p <- (1 + sum(permutation_moran >= observed_moran)) /
    (n_permutations + 1)
  moran_result <- data.frame(
    sample_id = sample_id,
    spots = length(viral_umi),
    knn_k = k_neighbors,
    undirected_edges = nrow(edges),
    library_depth_strata = length(index_groups),
    permutations = n_permutations,
    observed_moran_i_log1p = observed_moran,
    null_mean_moran_i = mean(permutation_moran),
    null_sd_moran_i = sd(permutation_moran),
    moran_z = (observed_moran - mean(permutation_moran)) / sd(permutation_moran),
    permutation_p_upper = moran_p,
    stringsAsFactors = FALSE
  )

  rm(object, counts)
  invisible(gc())
  list(threshold = threshold_results, moran = moran_result)
}

all_results <- lapply(
  seq_len(nrow(samples)),
  function(index) {
    analyse_sample(samples$sample_id[[index]], samples$input_file[[index]], index)
  }
)

threshold_results <- do.call(rbind, lapply(all_results, `[[`, "threshold"))
moran_results <- do.call(rbind, lapply(all_results, `[[`, "moran"))
moran_results$permutation_p_Bonferroni_2_samples <- pmin(
  1, 2 * moran_results$permutation_p_upper
)
threshold_results$permutation_p_maxT_then_Bonferroni_2_samples <- pmin(
  1, 2 * threshold_results$permutation_p_maxT_FWER
)
threshold_results$significant_clustered_experimentwise_0.05 <-
  threshold_results$join_count_z > 0 &
  threshold_results$permutation_p_maxT_then_Bonferroni_2_samples <= 0.05

write.table(
  threshold_results,
  file.path(output_dir, "Day4_PR8_spatial_randomness_threshold_scan.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  moran_results,
  file.path(output_dir, "Day4_PR8_spatial_randomness_MoranI.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

onset_rows <- do.call(
  rbind,
  lapply(samples$sample_id, function(sample_id) {
    candidate <- threshold_results[
      threshold_results$sample_id == sample_id &
        threshold_results$significant_clustered_experimentwise_0.05,
      , drop = FALSE
    ]
    if (nrow(candidate) == 0L) {
      return(data.frame(
        sample_id = sample_id,
        first_significant_threshold_raw_umi = NA_real_,
        first_significant_threshold_log1p = NA_real_,
        positive_spots_at_threshold = NA_integer_,
        positive_percent_at_threshold = NA_real_,
        experimentwise_p = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    candidate <- candidate[order(candidate$threshold_raw_umi), , drop = FALSE]
    data.frame(
      sample_id = sample_id,
      first_significant_threshold_raw_umi = candidate$threshold_raw_umi[[1L]],
      first_significant_threshold_log1p = candidate$threshold_log1p[[1L]],
      positive_spots_at_threshold = candidate$positive_spots[[1L]],
      positive_percent_at_threshold = candidate$positive_percent[[1L]],
      experimentwise_p =
        candidate$permutation_p_maxT_then_Bonferroni_2_samples[[1L]],
      stringsAsFactors = FALSE
    )
  })
)
write.table(
  onset_rows,
  file.path(output_dir, "Day4_PR8_first_nonrandom_threshold.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

threshold_plot <- ggplot(
  threshold_results,
  aes(
    x = threshold_raw_umi,
    y = -log10(permutation_p_maxT_then_Bonferroni_2_samples),
    colour = sample_id
  )
) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey35") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  scale_x_log10() +
  scale_colour_manual(values = c("2_2" = "#2166AC", "2_3" = "#67A9CF")) +
  labs(
    title = "Day 4 PR8-positive spots: spatial clustering across UMI thresholds",
    subtitle = paste0(
      k_neighbors, "-NN join-count test; ", n_permutations,
      " permutations within total-RNA-depth strata"
    ),
    x = "PR8-positive threshold (summed raw UMI per spot; log10 scale)",
    y = expression(-log[10](experimentwise~P)),
    colour = "Sample",
    caption = paste0(
      "Dashed line: P = 0.05. P values use maxT across thresholds, then ",
      "Bonferroni correction across two samples."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0)
  )

ggsave(
  file.path(output_dir, "Day4_PR8_spatial_randomness_threshold_scan.png"),
  threshold_plot, width = 9, height = 6.5, units = "in", dpi = 300, bg = "white"
)

report_lines <- c(
  "Day 4 PR8 spatial randomness test",
  "",
  "Samples: 2_2 and 2_3, tested separately.",
  paste0("PR8 expression: summed raw RNA UMI across detected PR8 features; thresholds use PR8_total_umi >= t."),
  paste0("Spatial graph: symmetrized ", k_neighbors, "-nearest-neighbour graph."),
  paste0(
    "Null: PR8 expression values are exchangeable among spots within each of ",
    library_strata_count, " rank-based total-RNA UMI strata."
  ),
  paste0("Permutations: ", n_permutations, "; random seed: ", random_seed, "."),
  paste0(
    "Continuous omnibus test: upper-tailed Moran's I on log1p(PR8 total UMI). ",
    "Two-sample Bonferroni P values are reported."
  ),
  paste0(
    "Threshold test: upper-tailed positive-positive join count. Candidate thresholds ",
    "must leave at least ", minimum_group_size, " positive and negative spots."
  ),
  paste0(
    "Threshold multiplicity: permutation maxT family-wise correction within each sample, ",
    "then Bonferroni correction across the two samples."
  ),
  "Interpretation: the first significant threshold is an operational detection boundary, not a biological infection cutoff or a causal transition.",
  "Fisher's exact test is not used because it ignores spatial adjacency; it remains suitable for a prespecified cell-type-by-positive 2x2 enrichment test."
)
writeLines(
  report_lines,
  file.path(output_dir, "README_Day4_PR8_spatial_randomness.txt")
)

message("Finished Day 4 PR8 spatial randomness tests.")
