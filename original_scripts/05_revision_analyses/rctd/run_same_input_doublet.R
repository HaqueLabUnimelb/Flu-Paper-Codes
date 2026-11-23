library(spacexr)

input_file <- paste0(
  "/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/",
  "RCTD/Day4mouseRCTDoutmin30_multi.rds"
)
output_file <- paste0(
  "/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/",
  "RCTD/Day4mouseRCTDoutmin30_doublet_same_input.rds"
)

message("R: ", R.version.string)
message("spacexr: ", as.character(packageVersion("spacexr")))
message("Loading same-input full object: ", input_file)

rctd_doublet <- readRDS(input_file)
if (!inherits(rctd_doublet, "RCTD.replicates")) {
  stop("Expected RCTD.replicates, got: ", paste(class(rctd_doublet), collapse = ", "))
}

message("Re-running identical RCTD inputs with doublet_mode='doublet'")
rctd_doublet <- run.RCTD.replicates(rctd_doublet, doublet_mode = "doublet")
saveRDS(rctd_doublet, output_file)
message("Saved same-input doublet result: ", output_file)
