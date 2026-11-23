#!/usr/bin/env Rscript
# Purpose: test whether two gene lists overlap more than expected by chance.
# Inputs: a background gene list and two gene lists (text files, one gene per line, first column is read).
# Output: the test printed to stdout, and a one-row tab-separated summary when an output file is given.
#
# Usage:
#   Rscript FishersExactTest.R BACKGROUND LIST_1 LIST_2 [OUTPUT_TSV]
# Paths can instead be supplied with FLU_FISHER_UNIVERSE, FLU_FISHER_GENE_LIST_1, FLU_FISHER_GENE_LIST_2
# and FLU_FISHER_OUTPUT.
#
# The calculation is that of original_scripts/05_revision_analyses/FishersExactTest.R: the number of
# background genes is the universe size N, the two lists are counted as read (a, b), their intersection is the
# overlap k, and a one-sided Fisher's exact test is run on matrix(c(k, a - k, b - k, N - (a + b - k)), nrow = 2).
# Genes of a tested list that are absent from the background file are counted as in the original; their number
# is reported in the summary.

args <- commandArgs(trailingOnly = TRUE)

path_argument <- function(position, env_name) {
  if (length(args) >= position && nzchar(args[[position]])) {
    return(args[[position]])
  }
  Sys.getenv(env_name, unset = "")
}

universe_file <- path_argument(1L, "FLU_FISHER_UNIVERSE")
list_1_file <- path_argument(2L, "FLU_FISHER_GENE_LIST_1")
list_2_file <- path_argument(3L, "FLU_FISHER_GENE_LIST_2")
output_file <- path_argument(4L, "FLU_FISHER_OUTPUT")

input_files <- c(universe_file, list_1_file, list_2_file)
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

total_gene_universe <- read.table(file = universe_file, stringsAsFactors = FALSE)[, 1]
gene_list_1 <- read.table(file = list_1_file, stringsAsFactors = FALSE)[, 1]
gene_list_2 <- read.table(file = list_2_file, stringsAsFactors = FALSE)[, 1]

# Inputs
A <- gene_list_1
B <- gene_list_2
N <- length(total_gene_universe)

k <- length(intersect(A, B))
a <- length(A)
b <- length(B)

# Build contingency table
mat <- matrix(c(
  k,
  a - k,
  b - k,
  N - (a + b - k)
), nrow = 2)

# Fisher's test
test <- fisher.test(mat, alternative = "greater")
print(test)

summary <- data.frame(
  universe_size = N,
  gene_list_1_size = a,
  gene_list_2_size = b,
  overlap_size = k,
  gene_list_1_outside_background = length(setdiff(A, total_gene_universe)),
  gene_list_2_outside_background = length(setdiff(B, total_gene_universe)),
  alternative = "greater",
  odds_ratio = unname(test$estimate),
  conf_low = test$conf.int[1],
  p_value = test$p.value,
  stringsAsFactors = FALSE
)
print(summary, digits = 17)

if (nzchar(output_file)) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.table(format(summary, digits = 17), output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved Fisher's exact-test summary: ", output_file)
}
