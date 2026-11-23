# RCTD revision analyses

Does the day-4 RCTD result depend on the RCTD mode, and where is each cell type?

## Order

1. Fit day 4 in full mode: `scripts/03_rctd/mouseDay4RCTD_full_v2.sh`.
2. `run_same_input_doublet.sh`: fits the same object again with `doublet_mode = "doublet"` (same beads,
   reference and normalisation).
3. `compare_full_vs_doublet.sh`: compares the two results.
4. `plot_individual_cell_types.sh`: one map per cell type and replicate (independent of 2 and 3).
5. `legacy/legacy_shared_spot_paired_ttest.sh`: compares the full-mode result with the earlier doublet-mode
   result `Day4mouseRCTDoutmin30_editv1`.

Submit from the root of the repository after `mkdir -p logs`, for example
`sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/05_revision_analyses/rctd/run_same_input_doublet.sh`.
The wrappers request one CPU, 16 to 32 GB (8 GB for the legacy comparison) and 31 days on `long`; override
these on the command line. The doublet-mode fit takes about 4 hours on one CPU and about 3 GB of memory.

## Configuration

- `FLU_RCTD_DAY4_FULL_RDS`: day-4 full-mode result (`RCTD.replicates` object);
- `FLU_RCTD_DAY4_DOUBLET_RDS`: doublet-mode result of step 2;
- `FLU_RCTD_DAY4_REPLICATE_NAMES`: replicate labels in object order (`puck1_lungST_2_2,puck2_lungST_2_3`);
- `FLU_RCTD_MODE_COMPARISON_OUTPUT_DIR`, `FLU_RCTD_CELL_PLOT_OUTPUT_DIR`: output folders;
- `FLU_RCTD_LEGACY_DOUBLET_RDS`, `FLU_RCTD_LEGACY_COMPARISON_OUTPUT_DIR`: for step 5.

Defaults are folders under `outputs/` in the repository.

## What the comparison computes

Bead classes of the doublet-mode fit; for the beads it does not reject: agreement of the first cell type,
the share of the full-mode weight in the two cell types of the doublet call, per-cell-type Pearson and
Spearman correlations, mean absolute and root mean square differences, total variation distance and
Jensen-Shannon divergence per bead, and mean compositions. Six summary plots, the tables and
`reviewer_summary.md` are written. Neighbouring beads are not independent, so these numbers describe agreement
between the two modes; they are not tests between mice.

## Comparison with the earlier result

The earlier doublet-mode result was fitted with `UMI_min = 75` and holds 16,292 and 14,421 beads, all of them
also in the full-mode result. `legacy_shared_spot_paired_ttest.R` compares the two results on the beads that the
doublet fit does not reject (9,983 and 8,977 beads), with a paired t-test per cell type. Because input and settings differ, a difference
cannot be assigned to the mode alone; step 3 compares the modes on the same input.
