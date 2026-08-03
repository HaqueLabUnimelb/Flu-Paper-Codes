# Analysis workflow

## 1. Curio Seeker preprocessing

`scripts/01_preprocessing/lung_curioseeker_qc.sh` launches the external Curio
Seeker Nextflow workflow with a samplesheet, output directory, work directory,
reference resource, Singularity profile, and Nextflow configuration. All paths
come from `config/paths.env`. Curio Seeker itself is not included.

## 2. Spatial feature plots

`scripts/02_feature_plots/FeaturePlots_Flu.R` reads four QCed Seurat objects,
normalizes each object with Seurat `NormalizeData()`, and generates:

- spatial `nFeature_RNA` QC plots;
- plots for Tnf, Ifng, Il6, Ifna1, Ifnb1, Cxcl1, Il10, and Il1b;
- plots for Igkc, Ighm, and Ighg2c.

The original decisions about which plots are saved are retained. Some feature
expressions are displayed but not saved by the script.

## 3. RCTD

Day 4 and day 10 are independent multi-replicate RCTD analyses. Both use two
spatial Seurat objects, a time-point-filtered single-cell reference, a minimum
UMI threshold of 30, one RCTD core, and full doublet mode. The recorded
`CELL_MIN_INSTANCE` values are 5 for day 4 and 2 for day 10.

The supplied SLURM wrappers request one CPU and 300 GB RAM. They enforce R 4.3.x,
Seurat 5.1.0, and spacexr 2.2.0.

## 4. RCTD visualization

`scripts/03_rctd/RCTD_flu_vis.R` normalizes each cell-type weight column to
the 0-to-1 range and overlays T cells, B cells, and epithelial cells for the
first replicate. It prints the plot and saves it only if
`FLU_RCTD_VIS_OUTPUT` is set.

## 5. Viral-reference mapping

The HKx31 and SARS-CoV-2 scripts launch the external Curio Seeker Nextflow
workflow with reference-specific samplesheets. These workflows are independent
of RCTD and can be submitted when their samplesheets and references are
available.

No full workflow was rerun during repository preparation.
