library(Seurat)
library(spacexr)

message("R: ", R.version.string)
message("Seurat: ", as.character(packageVersion("Seurat")))
message("spacexr: ", as.character(packageVersion("spacexr")))

# Seurat v5 can store counts in one or more assay layers. Return one sparse
# matrix while retaining cell names so annotations can be aligned explicitly.
get_counts <- function(object, assay = "RNA") {
  if (!assay %in% Assays(object)) {
    stop("Assay was not found: ", assay, ". Available assays: ",
         paste(Assays(object), collapse = ", "))
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

ref <- readRDS(
  "/data/projects/punim2394/lung_curio_mouse/lung_single_cell_ref/references/PR8_ref.rds"
)
ref <- ref[, ref$timepoint_ref == "10dpi"]
Idents(ref) <- ref$Annotation_relabel
ref$Annotation <- ref@active.ident
print(table(ref@active.ident))

puck_locs <- c(
  "/data/projects/punim2394/lung_curio_mouse_pr8/processed_data/lungST_4_2_QCed_edit.rds",
  "/data/projects/punim2394/lung_curio_mouse_pr8/processed_data/lungST_4_3_QCed_edit.rds"
)
rctd_save_loc <- paste0(
  "/data/projects/punim2394/lung_curio_mouse_pr8/",
  "OlliesFluRevisionCode/RCTD_day10"
)

pucks <- lapply(puck_locs, readRDS)
for (i in seq_along(pucks)) {
  counts <- get_counts(pucks[[i]], assay = "RNA")
  coords <- as.data.frame(Embeddings(pucks[[i]], reduction = "SPATIAL"))
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

reference <- Reference(
  counts = ref_count,
  cell_types = cell_types
)

if (identical(Sys.getenv("RCTD_PREFLIGHT_ONLY"), "1")) {
  message("Preflight reference: ", nrow(ref_count), " genes x ", ncol(ref_count), " cells")
  message("Preflight puck spots: ", paste(vapply(pucks, function(x) ncol(x@counts), integer(1)), collapse = ", "))
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
