# Original scripts

The files of the original analysis as they are on the Spartan cluster (University of Melbourne), copied byte
for byte: scripts, SLURM wrappers, samplesheets, configurations and run logs. Some scripts were saved again after
their outputs were written and some wrappers were not used (the script was run directly); both are noted in
[`../docs/PROVENANCE.md`](../docs/PROVENANCE.md).
[`MANIFEST.tsv`](MANIFEST.tsv) gives, for each file, its location on the cluster, owner, date, size and md5.
In six SLURM scripts the e-mail address of the `--mail-user` line is replaced by `<email>`; the manifest marks
them and gives the md5 of the source and of the copy.

| Folder | Content |
|---|---|
| `01_preprocessing/` | Curio Seeker v2.0.0 run that produced the Seurat objects: launcher, samplesheet, Nextflow configuration, log, and the script that builds the mouse + PR8 reference |
| `02_feature_plots/` | `FeaturePlots_Flu.R` as saved on 1 July 2026 (42 feature plots) and on 18 September 2026 (violin plots, second set of feature plots) |
| `03_rctd/` | RCTD of day 4 and day 10 (R scripts, SLURM wrappers, logs), the script that builds the conda environment with its package lists, and `RCTD_flu_vis.R` |
| `04_viral_mapping/` | Curio Seeker v3.1.0 runs with the mouse + X-31 and mouse + SARS-CoV-2 references: launchers, samplesheets, logs, Nextflow configuration, reference builds, and two R Markdown notebooks that read the results |
| `05_revision_analyses/` | Fisher's exact test of gene-list overlap; RCTD cell-type maps, doublet-mode fit of the same input and its comparison with full mode, comparison with the earlier doublet-mode result |
| `06_viral_expression/` | PR8 expression: maps, PR8-positive beads per cell type, Moran's I, distribution of counts, threshold series |

Runs, dates, software versions and references: [`../docs/PROVENANCE.md`](../docs/PROVENANCE.md). Versions with
configurable paths: [`../scripts/`](../scripts); their differences:
[`../docs/CHANGES_FROM_ORIGINAL.md`](../docs/CHANGES_FROM_ORIGINAL.md).
