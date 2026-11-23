# Analysis workflow

Versions, SLURM jobs and references of the original runs: [`PROVENANCE.md`](PROVENANCE.md).

## 1. Curio Seeker preprocessing

`scripts/01_preprocessing/lung_curioseeker_qc.sh` starts the Curio Seeker Nextflow workflow with a samplesheet,
an output folder, a work folder, the reference folder (`--igenomes_base`), the `singularity` profile and a
Nextflow configuration, all from `config/paths.env`. The original run used Curio Seeker v2.0.0 and a reference
of the mouse genome (mm10) joined with the PR8 genome, so that the PR8 genes are counted with the mouse genes.
Curio Seeker itself is not part of this repository.

## 2. Feature plots

`scripts/02_feature_plots/FeaturePlots_Flu.R` reads the four `lungST_*_QCed.rds` objects, normalises each with
Seurat `NormalizeData()` and draws on the tissue (`SPATIAL` reduction):

- the number of detected genes per bead (`nFeature_RNA`, log10 colour scale; files `mouse*_umis.jpeg`);
- Tnf, Ifng, Il6, Ifnb1, Cxcl1, Il10 and Il1b;
- Igkc, Ighm and Ighg2c.

Plots are saved as 12 cm x 12 cm JPEG files (magma colour scale). The second version of the script
(`original_scripts/02_feature_plots/FeaturePlots_Flu.R`) adds a violin plot per sample for the seven cytokine
genes and a second set of feature plots in which beads without counts are grey.

## 3. RCTD

Day 4 and day 10 are fitted separately, each with its two samples as replicates
(`create.RCTD.replicates`, `run.RCTD.replicates`): single-cell reference filtered to the matching time point,
`UMI_min = 30`, `max_cores = 1`, `doublet_mode = "full"`. `CELL_MIN_INSTANCE` is 5 for day 4 and 2 for day 10.
The wrappers check R 4.3, Seurat 5.1.0 and spacexr 2.2.0 before the fit.

`scripts/03_rctd/RCTD_flu_vis.R` scales every cell-type weight to 0-1 and overlays T cells, B cells and
epithelial cells of the first replicate.

## 4. Mapping to other viral references

`scripts/04_viral_mapping/` starts Curio Seeker with a reference of the mouse genome joined with the genome of
influenza A/X-31 (H3N2) or of SARS-CoV-2. The original runs used Curio Seeker v3.1.0.

## 5. Revision analyses

- `scripts/05_revision_analyses/FishersExactTest.R`: one-sided Fisher's exact test of the overlap of two gene
  lists, with the number of genes of a background list as the universe.
- `rctd/plot_individual_cell_types.R`: one map per cell type and day-4 replicate.
- `rctd/run_same_input_doublet.R` fits the day-4 result again with `doublet_mode = "doublet"` on the same
  beads and reference profiles; `rctd/compare_full_vs_doublet.R` compares the two results: bead classes,
  agreement of the first and second cell type, per-cell-type correlations and differences, total variation and
  Jensen-Shannon divergence per bead, mean compositions.
- `rctd/legacy/legacy_shared_spot_paired_ttest.R`: paired t-tests between the full-mode result and an earlier
  doublet-mode result, on the beads the two share.

## 6. PR8 expression

The scripts of `original_scripts/06_viral_expression/` sum the raw counts of the PR8 genes per bead and
describe this signal: maps, PR8-positive beads per cell type (the cell types stored in the annotated objects,
from an earlier doublet-mode RCTD run), spatial clustering (Moran's I and a join-count test), the distribution
of the counts and the dependence on the count threshold.
