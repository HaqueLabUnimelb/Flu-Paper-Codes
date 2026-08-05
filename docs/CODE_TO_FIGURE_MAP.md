# Code-to-figure map

No manuscript figure numbers were present in the supplied scripts, and no
figure-to-panel mapping was provided. The table records direct, code-supported
output relationships plus the authors' section assignment for the two viral
transcript PNG files whose plotting source script was not supplied.

| Code section | Generated output | Repository status | Manuscript mapping |
|---|---|---|---|
| `FeaturePlots_Flu.R`: spatial QC | `mouse1_umis.jpeg`, `mouse2_umis.jpeg`, `mouse4_umis.jpeg`, `mouse5_umis.jpeg` | 4 files included under `figures/generated/feature_plots/` | Unresolved |
| `FeaturePlots_Flu.R`: cytokine feature plots | Explicitly saved `mouse*_GENE.jpeg` cytokine plots | 27 files included under `figures/generated/feature_plots/` | Unresolved |
| `FeaturePlots_Flu.R`: B-cell-associated features | Explicitly saved `Igkc`, `Ighm`, and selected `Ighg2c` plots | 11 files included under `figures/generated/feature_plots/` | Unresolved |
| `RCTD_flu_vis.R`: cell-type overlay | Interactive plot; optional configured output | No existing, explicitly attributable file was found | Unresolved |
| Revision RCTD individual cell-type plots | One PNG per cell type and day-4 replicate plus a CSV summary | Code included; generated outputs not archived | Reviewer follow-up; mapping unresolved |
| Revision RCTD same-input mode comparison | Six PNG summaries, CSV metrics and `reviewer_summary.md` | Code included; generated outputs not archived | Reviewer follow-up; mapping unresolved |
| Revision Fisher gene-list overlap | TSV statistics only; no direct figure | Code included; result not archived | No direct figure output |
| Plotting source script not supplied: PR8 viral transcript mapping | `PR8_viral_gene_expression_2_2.png` and `PR8_viral_gene_expression_2_3.png` | 2 author-identified files included under `figures/generated/viral_transcript_mapping/` | Unresolved |

Curio Seeker, RCTD model-fitting, and viral-mapping launchers do not directly
open a graphics device or save a plot. External pipeline diagnostics were not
treated as direct outputs of the archived launch scripts.

All 44 included files were confirmed by exact filename, format, dimensions,
file size, and checksum. Direct script-to-output provenance is confirmed for the
42 feature-plot JPEG files. The two PR8 PNG files are assigned to the viral
transcript mapping section based on author-provided information because their
plotting source script was not supplied. No manuscript panel assignment is
inferred.
