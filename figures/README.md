# Archived figures

This directory contains 44 existing image outputs copied without modification.
The 42 JPEG feature-plot outputs have filenames and original output directories
that match explicit `ggsave()` calls in
`scripts/02_feature_plots/FeaturePlots_Flu.R`:

- 4 spatial `nFeature_RNA` QC plots;
- 27 explicitly saved cytokine-gene feature plots;
- 11 explicitly saved B-cell-associated gene feature plots.

The remaining two files are 2400 × 2100 RGB PNG outputs identified by the
authors as PR8 viral transcript mapping figures. The plotting source script was
not supplied, so their provenance is author-supplied rather than code-verified.

The JPEG files are stored in `generated/feature_plots/`, and the PNG files are
stored in `generated/viral_transcript_mapping/`. Together the 44 files occupy
28,387,196 bytes; the largest file is 3,172,863 bytes. Individual sizes and
SHA-256 checksums are recorded in `docs/FIGURE_MANIFEST.tsv` and
`figures/SHA256SUMS`.

No manuscript panel numbers were present in the scripts, so manuscript mappings
remain unresolved. The images were archived based on technical provenance and
the authors' stated section assignment only; no biological interpretation was
performed.
