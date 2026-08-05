#!/usr/bin/env Rscript
# Purpose: test whether two gene lists overlap more than expected by chance.
# Inputs: a background gene universe and two one-column gene-list text files.
# Output: a tab-separated summary written to FLU_FISHER_OUTPUT (and stdout).
#
# Usage:
#   Rscript FishersExactTest.R BACKGROUND LIST_1 LIST_2 [OUTPUT_TSV]
# Paths can instead be supplied with FLU_FISHER_UNIVERSE,
# FLU_FISHER_GENE_LIST_1, FLU_FISHER_GENE_LIST_2 and FLU_FISHER_OUTPUT.

args <- commandArgs(trailingOnly = TRUE)

path_argument <- function(position, env_name, default = "") {
  if (length(args) >= position && nzchar(args[[position]])) {
    return(args[[position]])
  }
  Sys.getenv(env_name, unset = default)
}

universe_file <- path_argument(1L, "FLU_FISHER_UNIVERSE")
list_1_file <- path_argument(2L, "FLU_FISHER_GENE_LIST_1")
list_2_file <- path_argument(3L, "FLU_FISHER_GENE_LIST_2")
output_file <- path_argument(4L, "FLU_FISHER_OUTPUT")

input_files <- c(universe_file, list_1_file, list_2_file)
input_labels <- c("background universe", "gene list 1", "gene list 2")
if (any(!nzchar(input_files))) {
  stop(
    "Supply BACKGROUND, LIST_1 and LIST_2 as command-line arguments or ",
    "set the corresponding FLU_FISHER_* environment variables."
  )
}
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop("Input file(s) not found: ", paste(missing_inputs, collapse = ", "))
}

read_gene_list <- function(path, label) {
  table <- read.table(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
  if (ncol(table) < 1L) {
    stop("No gene column was found in ", label, ": ", path)
  }
  genes <- trimws(as.character(table[[1L]]))
  genes <- unique(genes[nzchar(genes) & !is.na(genes)])
  if (length(genes) == 0L) {
    stop("No gene identifiers were found in ", label, ": ", path)
  }
  genes
}

gene_universe <- read_gene_list(universe_file, input_labels[[1L]])
gene_list_1 <- read_gene_list(list_1_file, input_labels[[2L]])
gene_list_2 <- read_gene_list(list_2_file, input_labels[[3L]])

outside_universe <- list(
  gene_list_1 = setdiff(gene_list_1, gene_universe),
  gene_list_2 = setdiff(gene_list_2, gene_universe)
)
outside_counts <- lengths(outside_universe)
if (any(outside_counts > 0L)) {
  details <- paste(
    names(outside_counts)[outside_counts > 0L],
    outside_counts[outside_counts > 0L],
    sep = "=",
    collapse = ", "
  )
  stop("Every tested gene must be in the background universe; outside counts: ", details)
}

overlap_size <- length(intersect(gene_list_1, gene_list_2))
list_1_size <- length(gene_list_1)
list_2_size <- length(gene_list_2)
universe_size <- length(gene_universe)

# Rows represent membership in list 1; columns represent membership in list 2.
contingency_table <- matrix(
  c(
    overlap_size,
    list_1_size - overlap_size,
    list_2_size - overlap_size,
    universe_size - list_1_size - list_2_size + overlap_size
  ),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(
    list_1 = c("in", "not_in"),
    list_2 = c("in", "not_in")
  )
)
if (any(contingency_table < 0L)) {
  stop("The gene lists and background universe produced an invalid contingency table.")
}

test <- fisher.test(contingency_table, alternative = "greater")
summary <- data.frame(
  universe_size = universe_size,
  gene_list_1_size = list_1_size,
  gene_list_2_size = list_2_size,
  overlap_size = overlap_size,
  alternative = "greater",
  odds_ratio = unname(test$estimate),
  p_value = test$p.value,
  stringsAsFactors = FALSE
)

write.table(summary, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
if (nzchar(output_file)) {
  output_dir <- dirname(output_file)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(summary, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved Fisher's exact-test summary: ", output_file)
}
