#!/usr/bin/env Rscript
# Purpose: run multi-replicate RCTD for the two day-10 spatial samples.
# Inputs: FLU_REFERENCE_RDS, FLU_DAY10_PUCK1_RDS and FLU_DAY10_PUCK2_RDS.
# Output: Day10mouseRCTDoutmin30_multi.rds in FLU_RCTD_OUTPUT_DIR.
# Scientific parameters from the archived analysis are retained below.

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
})

message("R: ", R.version.string)
message("Seurat: ", as.character(packageVersion("Seurat")))
message("spacexr: ", as.character(packageVersion("spacexr")))

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Unable to determine the script location. Run this file with Rscript.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", file_argument)))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) {
    stop("Required path is empty: ", name)
  }
  value
}

reference_file <- env_path(
  "FLU_REFERENCE_RDS",
  file.path(repo_dir, "data", "processed", "PR8_ref.rds")
)
puck_locs <- c(
  env_path(
    "FLU_DAY10_PUCK1_RDS",
    file.path(repo_dir, "data", "processed", "lungST_4_2_QCed_edit.rds")
  ),
  env_path(
    "FLU_DAY10_PUCK2_RDS",
    file.path(repo_dir, "data", "processed", "lungST_4_3_QCed_edit.rds")
  )
)
rctd_save_loc <- env_path(
  "FLU_RCTD_OUTPUT_DIR",
  file.path(repo_dir, "outputs", "rctd")
)

missing_inputs <- c(reference_file, puck_locs)[!file.exists(c(reference_file, puck_locs))]
if (length(missing_inputs) > 0L) {
  stop("Required RDS input(s) not found: ", paste(missing_inputs, collapse = ", "))
}
dir.create(rctd_save_loc, recursive = TRUE, showWarnings = FALSE)

# Seurat v5 can store counts in one or more assay layers. Return one sparse
# matrix while retaining cell names so annotations can be aligned explicitly.
get_counts <- function(object, assay = "RNA") {
  if (!assay %in% Assays(object)) {
    stop(
      "Assay was not found: ", assay, ". Available assays: ",
      paste(Assays(object), collapse = ", ")
    )
  }

  assay_object <- object[[assay]]
  if (!inherits(assay_object, "Assay5")) {
    return(as(LayerData(assay_object, layer = "counts"), "dgCMatrix"))
  }

  count_layers <- Layers(assay_object, search = "^counts")
  if (length(count_layers) == 0L) {
    stop("No counts layer was found in assay: ", assay)
  }

  count_matrices <- lapply(count_layers, function(layer) {
    as(LayerData(assay_object, layer = layer), "dgCMatrix")
  })
  common_genes <- Reduce(intersect, lapply(count_matrices, rownames))
  if (length(common_genes) == 0L) {
    stop("The counts layers do not share any genes.")
  }
  count_matrices <- lapply(
    count_matrices,
    function(x) x[common_genes, , drop = FALSE]
  )
  counts <- do.call(cbind, count_matrices)

  if (anyDuplicated(colnames(counts))) {
    stop("Duplicated cell/spot names were found after joining counts layers.")
  }
  counts
}

ref <- readRDS(reference_file)
ref <- ref[, ref$timepoint_ref == "10dpi"]
Idents(ref) <- ref$Annotation_relabel
ref$Annotation <- ref@active.ident
print(table(ref@active.ident))

pucks <- lapply(puck_locs, readRDS)
for (i in seq_along(pucks)) {
  counts <- get_counts(pucks[[i]], assay = "RNA")
  coords <- as.data.frame(Embeddings(pucks[[i]], reduction = "SPATIAL"))
  if (ncol(coords) != 2L) {
    stop("The SPATIAL reduction must contain exactly two coordinate columns.")
  }
  colnames(coords) <- c("x", "y")

  missing_coords <- setdiff(colnames(counts), rownames(coords))
  if (length(missing_coords) > 0L) {
    stop("Missing SPATIAL coordinates for ", length(missing_coords), " spots.")
  }
  coords <- coords[colnames(counts), , drop = FALSE]
  nUMI <- Matrix::colSums(counts)

  pucks[[i]] <- new(
    "SpatialRNA",
    counts = counts,
    coords = coords,
    nUMI = nUMI
  )
}

ref_count <- get_counts(ref, assay = "RNA")
cell_types <- ref@active.ident[colnames(ref_count)]
if (anyNA(cell_types)) {
  stop("Cell type annotations could not be aligned to all reference cells.")
}
cell_types <- factor(as.character(cell_types))
names(cell_types) <- colnames(ref_count)

reference <- Reference(counts = ref_count, cell_types = cell_types)

if (identical(Sys.getenv("RCTD_PREFLIGHT_ONLY"), "1")) {
  message(
    "Preflight reference: ", nrow(ref_count), " genes x ",
    ncol(ref_count), " cells"
  )
  message(
    "Preflight puck spots: ",
    paste(vapply(pucks, function(x) ncol(x@counts), integer(1)), collapse = ", ")
  )
  message("RCTD input compatibility preflight passed.")
  quit(save = "no", status = 0L)
}

rctd.obj <- create.RCTD.replicates(
  pucks,
  reference,
  max_cores = 1,
  CELL_MIN_INSTANCE = 2,
  replicate_names = c("puck1", "puck2"),
  UMI_min = 30
)

rm(pucks, ref, reference, ref_count, cell_types)
gc()

rctd.obj <- run.RCTD.replicates(rctd.obj, doublet_mode = "full")
output_file <- file.path(rctd_save_loc, "Day10mouseRCTDoutmin30_multi.rds")
saveRDS(rctd.obj, output_file)
message("Saved RCTD result to: ", output_file)
