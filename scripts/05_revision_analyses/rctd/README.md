# RCTD revision follow-up analyses

These scripts evaluate whether the day-4 RCTD result is sensitive to full versus
doublet mode and provide per-cell-type spatial maps.

## Recommended order

1. Generate the day-4 full-mode object with
   `scripts/03_rctd/mouseDay4RCTD_full_v2.sh`.
2. Submit `run_same_input_doublet.sh` to rerun that fitted object with
   `doublet_mode = "doublet"`.
3. Submit `compare_full_vs_doublet.sh` to calculate spot- and cell-type-level
   agreement metrics and reviewer-facing plots.
4. Submit `plot_individual_cell_types.sh` only when individual spatial
   proportion maps are required.

Create `logs/` before calling `sbatch`. The wrappers use the `long`
partition for 31 days, one CPU and 16--32 GB RAM. They intentionally omit
`--account`; supply an authorized account at submission time if required:

```bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT scripts/05_revision_analyses/rctd/run_same_input_doublet.sh
sbatch --account=YOUR_ACCOUNT scripts/05_revision_analyses/rctd/compare_full_vs_doublet.sh
sbatch --account=YOUR_ACCOUNT scripts/05_revision_analyses/rctd/plot_individual_cell_types.sh
```

## Configuration

The main variables are:

- `FLU_RCTD_DAY4_FULL_RDS`: day-4 full-mode `RCTD.replicates` object;
- `FLU_RCTD_DAY4_DOUBLET_RDS`: same-input doublet-mode result;
- `FLU_RCTD_DAY4_REPLICATE_NAMES`: comma-separated labels in object order;
- `FLU_RCTD_MODE_COMPARISON_OUTPUT_DIR`: comparison tables and plots;
- `FLU_RCTD_CELL_PLOT_OUTPUT_DIR`: per-cell-type maps.

Defaults are repository-relative under `outputs/`.

## Comparison outputs

The same-input comparison writes per-replicate and per-cell-type metrics,
spot-class counts, spot-level agreement, mean composition tables, six PNG
summaries and `reviewer_summary.md`. Doublet-mode reject spots are excluded
from proportion agreement metrics. The primary summaries are top-cell-type
concordance, full-mode mass in the top two cell types, correlations, mean
absolute error, total variation distance and Jensen-Shannon divergence.

Spot-level observations are spatially autocorrelated. These metrics describe
agreement between modes and are not a substitute for biological-replicate
inference.

## Legacy analysis

The [legacy directory](legacy/README.md) contains an earlier shared-spot paired
t-test analysis for provenance. It is not the recommended mode-sensitivity
analysis.
