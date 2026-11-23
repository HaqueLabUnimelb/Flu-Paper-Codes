# Provenance: data, runs and software

What produced the results, taken from the run logs, the SLURM accounting records and the files on the
Spartan cluster (University of Melbourne), as checked on 23 November 2026. Paths of the cluster appear as they
are in [`../original_scripts/`](../original_scripts) and its [`MANIFEST.tsv`](../original_scripts/MANIFEST.tsv).

## Samples and data

Lungs of four female C57BL/6J mice after intranasal infection with 500 PFU of influenza A/Puerto Rico/8/1934
(PR8, H1N1; PR8-OVA in the GEO record), profiled with Curio Seeker (10 µm beads). GEO series [GSE341918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341918)
(public since 15 September 2026).

| Sample in the scripts | Library | Day after infection | GEO sample | Bead barcode file | Object on GEO (RCTD input) | Beads x genes |
|---|---|---|---|---|---|---|
| mouse1 | lungST_2.2 | 4 | GSM9920041 | A0090_063 | `lungST_2_2_QCed_edit.rds` | 42,390 x 21,226 |
| mouse2 | lungST_2.3 | 4 | GSM9920042 | A0090_067 | `lungST_2_3_QCed_edit.rds` | 47,207 x 21,439 |
| mouse4 | lungST_4.2 | 10 | GSM9920043 | A0090_061 | `lungST_4_2_QCed_edit.rds` | 41,816 x 18,750 |
| mouse5 | lungST_4.3 | 10 | GSM9920044 | A0090_062 | `lungST_4_3_QCed_edit.rds` | 30,616 x 18,402 |

Three sets of Seurat objects are read by the scripts, all from the same preprocessing run:

| Objects | Beads (2.2 / 2.3 / 4.2 / 4.3) | Read by |
|---|---|---|
| `lungST_*_QCed_edit.rds` (on GEO) | 42,390 / 47,207 / 41,816 / 30,616 | RCTD |
| `lungST_*_QCed.rds` | 42,113 / 46,545 / 39,900 / 28,823 | feature plots and violin plots |
| `annotated_lungST_*_QCed.rds` (with the cell types of an earlier doublet-mode RCTD run) | 31,666 / 29,053 / 16,098 / 12,517 | PR8 expression analyses |

Each object has one `RNA` assay and a two-dimensional `SPATIAL` reduction. The PR8 genes are features of the
objects: `HA, M1, M2, NA, NEP, NP, NS1, PA, PB1, PB2` on day 4; nine and six of them on day 10.

## 1. Preprocessing

| Item | Value |
|---|---|
| Pipeline | Curio Seeker v2.0.0 (container `curio-seeker-singularity:v2.0.0`), Nextflow 23.10.1, Singularity |
| Run | SLURM job 1634899, 15 November 2024, 2 h 32 min; log `original_scripts/01_preprocessing/slurm-1634899.out` |
| Launcher, samplesheet, configuration | `original_scripts/01_preprocessing/` |
| Reference | genome key `mm10_opt_PR8`: mouse mm10 (`mouse_mm10_optimized_v1` FASTA and GTF) joined with the PR8 genome, RefSeq assembly GCF_000865725.1 (segments NC_002016.1 to NC_002023.1); STAR 2.6.1d index. Built with `generate_mouse_PR8_reference.sh` and Curio's `generate_seeker_reference.sh` |
| Output used | `lungST_{2.2,2.3,4.2,4.3}_seurat.rds` (written by the pipeline step `FORMATCLEANUP`). The later pipeline steps of this run (h5ad conversion, report) ended with an error and their outputs are not used |

## 2. Single-cell reference for RCTD

`PR8_ref.rds`: 17,826 cells from two public data sets of PR8-infected mouse lungs, matched to the deposited
files by cell barcode.

| Part of the reference | Cells | Source | `timepoint_ref` |
|---|---:|---|---|
| 3 days after infection (one sample, `Day3_young3`) | 5,777 | GEO GSE202325 | `4dpi` |
| 5 days, immune and epithelial samples | 2,911 | GEO GSE184384 | `4dpi` |
| 9 and 11 days, immune and epithelial samples | 9,138 | GEO GSE184384 | `10dpi` |

The day-4 fit uses the 8,688 cells with `timepoint_ref == "4dpi"` (12 cell types of `Annotation_relabel`), the
day-10 fit the 9,138 cells with `"10dpi"` (10 cell types).

## 3. RCTD

spacexr 2.2.0 (source: commit `205c9ac` of <https://github.com/dmcable/spacexr>), R 4.3.3, Seurat 5.1.0,
SeuratObject 5.2.0, in the conda environment built by `original_scripts/03_rctd/setup_spacexr2.2.sh`
(`environment_spacexr2.2.yml`; all conda packages in `conda_spacexr2.2_explicit.txt`; spacexr installed by SLURM
job 28239253, 28 July 2026, log `original_scripts/03_rctd/logs/rctd22_setup_28239253.log`). Both fits:
`create.RCTD.replicates(..., UMI_min = 30, max_cores = 1)` and `run.RCTD.replicates(doublet_mode = "full")`,
two replicates per time point.

| | Day 4 | Day 10 |
|---|---|---|
| SLURM job | 28239254, 28 July 2026 | 28289107, 29 July 2026 |
| Wrapper and script | `mouseDay4RCTD_full_v2.sh`, `RCTD_flu_d4_multi_spacexr2.2.R` | `mouseDay10RCTD_full_v3.sh`, `RCTD_flu_d10_multi_spacexr2.2.R` |
| Run time, peak memory | 38 min 31 s, 4.5 GiB (300 GB requested) | 15 min 24 s, 7.8 GiB (300 GB requested) |
| Reference | 8,688 cells, 12 cell types | 9,138 cells, 10 cell types |
| `CELL_MIN_INSTANCE` | 5 | 2 (the reference holds two endothelial cells at this time point) |
| Beads with weights (replicate 1, 2) | 31,644 and 29,295 | 15,362 and 11,560 |
| Result | `Day4mouseRCTDoutmin30_multi.rds` | `Day10mouseRCTDoutmin30_multi.rds` |
| Log | `original_scripts/03_rctd/logs/day4_rctd22_28239254.log` | `.../logs/day10_rctd_28289107.log` |

The day-4 R script was saved once more on 29 July 2026, after its run: the run wrote
`Day4mouseRCTDoutmin30_multi`, the saved script writes `Day4mouseRCTDoutmin30_multi.rds`. The settings stored in
the result object are the ones above.

`RCTD_flu_vis.R` overlays the weights of T cells, B cells and epithelial cells of the first replicate. It names
the epithelial type `Epithelial Cells`, the label of the earlier doublet-mode results; in the two full-mode
results above the label is `Epithelial cells`.

## 4. Feature plots and violin plots

| Version of `FeaturePlots_Flu.R` | Output | Written |
|---|---|---|
| saved 1 July 2026 (`original_scripts/02_feature_plots/FeaturePlots_Flu_2026-07-01.R`) | 42 JPEG files: `nFeature_RNA` per bead (`mouse*_umis.jpeg`), cytokine genes, B-cell genes; `figures/generated/feature_plots/` | 1 July 2026 |
| saved 18 September 2026 (`.../FeaturePlots_Flu.R`) | 7 violin plots (`figures/generated/violin_plots/`) and 27 feature plots with undetected beads in grey (`figures/generated/feature_plots_v2/`) | 9 and 17 August 2026 |

The scripts were run in interactive RStudio sessions. The 42 JPEG files of the first version are reproduced
byte for byte with R 4.3.1, Seurat 4.3.0.1, SeuratObject 4.1.3, ggplot2 3.5.0, tidyverse 2.0.0 and
mltools 0.3.5 (see [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md)). The second version was saved after its plots
were written; for `Ifnb1` of mouse4 it saves a JPEG file, while the folder of the plots holds a PDF file of that
plot.

## 5. Mapping to other viral references

The reads of the four lungs were also mapped to the mouse genome joined with the genome of influenza
A/X-31 (H3N2) or of SARS-CoV-2.

| Run | Samples | SLURM job | Reference |
|---|---|---|---|
| X-31, day 4 | lungST_2.2, 2.3 | 19877908, 18 December 2025 | `GRCm38_HKX31` |
| X-31, day 10 | lungST_4.2, 4.3 | 19887314, 19 December 2025 | `GRCm38_HKX31` |
| SARS-CoV-2 | all four | 19891965, 19 December 2025 | `GRCm38_covid19` |

Pipeline: Curio Seeker v3.1.0, Nextflow 23.04.2, Apptainer 1.4.4. References: mouse GRCm38 with the GENCODE
vM25 annotation, joined with the X-31 segments GenBank DQ874873.1 to DQ874880.1, or with SARS-CoV-2 MN908947.3
(Ensembl ASM985889v3, release 101); built with `reference.slurm` / `reference_covid.slurm`, Curio's
`generate_seeker_reference.sh` and STAR 2.6.1d. Launchers, samplesheets, logs and the two R Markdown notebooks
that read the results are in `original_scripts/04_viral_mapping/`.

## 6. PR8 expression

Scripts of `original_scripts/06_viral_expression/`, run on 3 to 11 August 2026 with the R environment of the
RCTD fits (R 4.3.3, Seurat 5.1.0, ggplot2 3.5.2; `sessionInfo.txt`) on the `annotated_lungST_*_QCed.rds`
objects: three through their wrappers `run_*.slurm` (spatial pattern, distribution, threshold series), the others
directly with `Rscript`. The cell types are those stored in the annotated objects (an earlier doublet-mode RCTD
run, 14 cell types). `assess_pr8_viral_expression.R` was saved again at 16:44 on 3 August, after its outputs
(14:50).

| Script | What it does | Output in `figures/generated/` |
|---|---|---|
| `assess_pr8_viral_expression.R` | per bead, sum of the raw counts of the PR8 genes; map of `log1p` of the sum, one colour scale for the four samples (upper limit: 99th percentile of the positive beads, 51 UMI) | `viral_transcript_mapping/PR8_viral_gene_expression_*.png`, `viral_expression/pr8_*` |
| `plot_day4_rctd_infected_spots.R`, `plot_day4_rctd_infected_bar_with_percentage_line.R` | PR8-positive beads per cell type, day 4 | `viral_expression/Day4_RCTD_celltype_infection_*`, `Day4_RCTD_celltype_infected_bar_with_percentage_line.png` |
| `test_day4_pr8_spatial_randomness.R` | spatial clustering of the PR8 signal, day 4: Moran's I of `log1p` of the counts, and a join-count test of PR8-positive beads at count thresholds (permutations within strata of total counts) | `viral_expression/Day4_PR8_*` |
| `plot_pr8_8segment_violin.R` | distribution of the PR8 counts per bead, four samples | `viral_expression/PR8_8segment_*` |
| `plot_day4_rctd_threshold_series.R`, `plot_threshold_celltype_lineplot.R` | share of PR8-positive beads per cell type at 14 count thresholds from 1 to 50 | `viral_expression/Day4_RCTD_celltype_raw_UMI_threshold_series.tsv`, `Day4_RCTD_infected_percentage_thresholds_*` |

Four of these scripts write a note `README_*.txt` with their definitions (in `figures/generated/viral_expression/`).

## 7. Revision analyses

| Analysis | Script | Run |
|---|---|---|
| Overlap of two gene lists | `FishersExactTest.R` | Background list `backgroundgenelist.txt` (21,214 genes); SpatialDE gene lists of mouse 4 (61 genes) and mouse 5 (48 genes), six genes of each not in the background list. Overlap 29 genes, odds ratio 1,024, one-sided p = 7.4e-65 |
| Map of every RCTD cell type, day 4 | `rctd/plot_individual_cell_types.R` | SLURM job 28261808, 29 July 2026: 24 maps |
| Full mode against doublet mode, same input | `rctd/run_same_input_doublet.R`, `rctd/compare_full_vs_doublet.R` | the day-4 result is fitted again with `doublet_mode = "doublet"` and compared bead by bead. SLURM job 28262524, 29 July 2026: the fit ran (3 h 46 min) and the job stopped when it wrote the result. Fit and comparison were run on 23 November 2026 ([`REPRODUCIBILITY.md`](REPRODUCIBILITY.md)) |
| Full mode against the earlier doublet-mode result | `rctd/legacy_shared_spot_paired_ttest.R` | SLURM job 28266050, 29 July 2026. The earlier result (doublet mode, `UMI_min = 75`) holds 16,292 and 14,421 beads, all of them also in the full-mode result; the paired t-tests use the 9,983 and 8,977 of them that the doublet fit does not reject |
