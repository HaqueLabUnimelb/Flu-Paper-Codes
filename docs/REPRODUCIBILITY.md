# Tests on the original data

On 23 November 2026 the R scripts of `scripts/` were run on the Spartan cluster on the original inputs, those
with a SLURM wrapper through it, and their outputs were compared with the outputs of the original runs.

| Script | Software | Result |
|---|---|---|
| `03_rctd/mouseDay4RCTD_full_v2.sh` (day 4) | R 4.3.3, Seurat 5.1.0, spacexr 2.2.0 | Weights of both replicates identical to the original result (31,644 x 12 and 29,295 x 12 beads x cell types; largest difference 0); differentially expressed gene counts, likelihoods and sigma equal to the original log. 40 min |
| `03_rctd/mouseDay10RCTD_full_v3.sh` (day 10) | as above | Weights identical (15,362 x 10 and 11,560 x 10; largest difference 0). 16 min |
| `05_revision_analyses/rctd/plot_individual_cell_types.sh` | as above | 24 maps identical byte for byte; summary table identical apart from its column of output paths |
| `05_revision_analyses/rctd/legacy/legacy_shared_spot_paired_ttest.sh` | as above | the four tables and the note identical byte for byte |
| `05_revision_analyses/rctd/run_same_input_doublet.sh` | as above | the doublet-mode fit of the day-4 input ran in 3 h 51 min with 3.1 GB and wrote its result (31,644 and 29,295 beads, the beads of the full-mode result) |
| `05_revision_analyses/rctd/compare_full_vs_doublet.sh` | as above | on that result, the five tables, six plots and the summary are identical byte for byte to those of the original script (paths changed) run on the same input. Beads not rejected by the doublet fit: 17,841 of 31,644 and 14,768 of 29,295; same first cell type in 74.6 % and 76.9 % of them; median Spearman correlation of a cell type's proportions between the modes 0.18 and 0.26 |
| `05_revision_analyses/FishersExactTest.R` | R 4.3.3 | the original script's numbers: universe 21,214 genes, lists of 61 and 48 genes, overlap 29, odds ratio 1024.2638360311064, p = 7.4440129213981662e-65 |
| `02_feature_plots/FeaturePlots_Flu.R` | R 4.3.1, Seurat 4.3.0.1, SeuratObject 4.1.3, ggplot2 3.5.0, patchwork 1.2.0, tidyverse 2.0.0, data.table 1.14.8, mltools 0.3.5, spacexr 2.2.1 | the 42 JPEG files of `figures/generated/feature_plots/` identical byte for byte |
| same | R 4.5.0, Seurat 5.5.1, ggplot2 4.0.3 | the same 42 files with the same layout and colour scale; points are drawn larger |
| `03_rctd/RCTD_flu_vis.R` | both software sets above | runs on the two full-mode results and on the two earlier doublet-mode results; the plot is written when `FLU_RCTD_VIS_OUTPUT` is set |

Further checks:

- The wrappers were submitted with `sbatch` from the root of the repository (cell-type maps, day-4 RCTD input
  check, day-10 RCTD) and from another directory with `FLU_REPO_DIR` exported (comparison with the earlier
  result). The day-4 input check gave a reference of 25,954 genes x 8,688 cells and 42,390 and 47,207 beads, as
  in the original run.
- The Curio Seeker launchers, which need the FASTQ files and the pipeline containers, were checked by their
  command line: with the paths of the original runs filled in, they give the same `nextflow run` options as the
  original launchers.
- Repeated full-mode RCTD runs give identical weights. spacexr draws a random subset of beads to estimate sigma;
  in both fitting scripts the subset of the reference object fixes the state of R's random number generator before
  that step. The doublet-mode refit has no such step, so repeated runs differ slightly (sigma of replicate 1: 0.84
  in the original job, 0.82 in the test run), and so do the numbers of the comparison.
- `tests/validate_repository.sh` passes all checks.

The test records (job logs, comparison tables) are kept on the lab's storage.
