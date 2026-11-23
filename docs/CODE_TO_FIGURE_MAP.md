# Scripts and their outputs

Every file of `figures/generated/` with the script that wrote it, its original location, date and checksum is
listed in [`FIGURE_MANIFEST.tsv`](FIGURE_MANIFEST.tsv). The files are the outputs of the original runs,
unchanged.

| Script (`original_scripts/`) | Output | In `figures/generated/` |
|---|---|---|
| `02_feature_plots/FeaturePlots_Flu_2026-07-01.R` | detected genes per bead (`mouse*_umis.jpeg`, 4 files); Tnf, Ifng, Il6, Ifnb1, Cxcl1, Il10, Il1b (27 files); Igkc, Ighm, Ighg2c (11 files) | `feature_plots/` (42) |
| `02_feature_plots/FeaturePlots_Flu.R` | violin plot per sample of Il1b, Il6, Cxcl1, Tnf, Ifng, Ifnb1, Il10 | `violin_plots/` (7) |
| `02_feature_plots/FeaturePlots_Flu.R` | the same seven genes on the tissue, beads without counts in grey | `feature_plots_v2/` (27) |
| `06_viral_expression/assess_pr8_viral_expression.R` | PR8 counts per bead on the tissue, four samples; summary tables | `viral_transcript_mapping/` (4); `viral_expression/pr8_*`, `README_results.txt` |
| `06_viral_expression/plot_day4_rctd_infected_spots.R` | PR8-positive beads per cell type, day 4 | `viral_expression/Day4_RCTD_celltype_infection_*`, `README_Day4_RCTD_celltype_infection.txt` |
| `06_viral_expression/plot_day4_rctd_infected_bar_with_percentage_line.R` | the same as a bar plot with the percentage | `viral_expression/Day4_RCTD_celltype_infected_bar_with_percentage_line.png` |
| `06_viral_expression/test_day4_pr8_spatial_randomness.R` | spatial clustering of the PR8 signal (Moran's I, join-count test) | `viral_expression/Day4_PR8_*`, `README_Day4_PR8_spatial_randomness.txt` |
| `06_viral_expression/plot_pr8_8segment_violin.R` | PR8 counts per bead, distribution | `viral_expression/PR8_8segment_*`, `README_PR8_8segment_violin.txt` |
| `06_viral_expression/plot_day4_rctd_threshold_series.R` | PR8-positive share per cell type at 14 count thresholds | `viral_expression/Day4_RCTD_celltype_raw_UMI_threshold_series.tsv` |
| `06_viral_expression/plot_threshold_celltype_lineplot.R` | the same at thresholds 1, 2, 5 and 10 as line plots | `viral_expression/Day4_RCTD_infected_percentage_thresholds_*` |

Scripts whose outputs are data objects or sets of maps:

| Script | Output |
|---|---|
| `01_preprocessing/lung_curioseeker_qc.sh` | Curio Seeker output folder with one Seurat object per sample |
| `03_rctd/RCTD_flu_d4_multi_spacexr2.2.R`, `RCTD_flu_d10_multi_spacexr2.2.R` | `Day4mouseRCTDoutmin30_multi.rds`, `Day10mouseRCTDoutmin30_multi.rds` |
| `03_rctd/RCTD_flu_vis.R` | overlay of the T-cell, B-cell and epithelial weights, drawn on screen; written for the earlier doublet-mode results (label `Epithelial Cells`). `scripts/03_rctd/RCTD_flu_vis.R` also reads the full-mode results (label `Epithelial cells`) and saves the plot |
| `04_viral_mapping/` launchers | Curio Seeker output folders of the X-31 and SARS-CoV-2 runs |
| `05_revision_analyses/FishersExactTest.R` | the test, printed |
| `05_revision_analyses/rctd/plot_individual_cell_types.R` | 24 maps (12 cell types, two day-4 replicates) and `cell_type_plot_summary.csv` |
| `05_revision_analyses/rctd/run_same_input_doublet.R`, `compare_full_vs_doublet.R` | doublet-mode result of the day-4 input; tables, six summary plots and `reviewer_summary.md` of the comparison (results of 23 November 2026 in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md)) |
| `05_revision_analyses/rctd/legacy_shared_spot_paired_ttest.R` | proportions of the compared beads in both results and the paired t-tests |
