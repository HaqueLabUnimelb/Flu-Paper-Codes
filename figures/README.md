# Output figures and tables

`generated/` holds 107 files written by the original runs, copied without change. Script, original location,
date, size and SHA-256 of each file: [`../docs/FIGURE_MANIFEST.tsv`](../docs/FIGURE_MANIFEST.tsv); checksums
also in [`SHA256SUMS`](SHA256SUMS) (`sha256sum -c figures/SHA256SUMS` from the root of the repository).

| Folder | Files | Content | Script |
|---|---:|---|---|
| `generated/feature_plots/` | 42 | detected genes per bead (`mouse*_umis.jpeg`), cytokine genes and B-cell genes on the tissue | `FeaturePlots_Flu_2026-07-01.R` |
| `generated/violin_plots/` | 7 | expression of seven cytokine genes per sample | `FeaturePlots_Flu.R` |
| `generated/feature_plots_v2/` | 27 | the seven cytokine genes on the tissue, beads without counts in grey | `FeaturePlots_Flu.R` |
| `generated/viral_transcript_mapping/` | 4 | PR8 counts per bead on the tissue, one map per sample | `assess_pr8_viral_expression.R` |
| `generated/viral_expression/` | 27 | PR8 signal per cell type, its spatial pattern and its distribution: plots, the tables behind them and the notes written by the scripts | scripts of `original_scripts/06_viral_expression/` |

Samples: mouse1 = lungST_2.2 and mouse2 = lungST_2.3 (day 4); mouse4 = lungST_4.2 and mouse5 = lungST_4.3
(day 10).
