# Code-to-figure map

No manuscript figure numbers were present in the supplied scripts, and no
figure-to-panel mapping was provided. The table therefore records only direct,
code-supported output relationships.

| Code section | Generated output | Repository status | Manuscript mapping |
|---|---|---|---|
| `FeaturePlots_Flu.R`: spatial QC | `mouse1_umis.jpeg`, `mouse2_umis.jpeg`, `mouse4_umis.jpeg`, `mouse5_umis.jpeg` | 4 files included under `figures/generated/feature_plots/` | Unresolved |
| `FeaturePlots_Flu.R`: cytokine feature plots | Explicitly saved `mouse*_GENE.jpeg` cytokine plots | 27 files included under `figures/generated/feature_plots/` | Unresolved |
| `FeaturePlots_Flu.R`: B-cell-associated features | Explicitly saved `Igkc`, `Ighm`, and selected `Ighg2c` plots | 11 files included under `figures/generated/feature_plots/` | Unresolved |
| `RCTD_flu_vis.R`: cell-type overlay | Interactive plot; optional configured output | No existing, explicitly attributable file was found | Unresolved |

Curio Seeker, RCTD model-fitting, and viral-mapping launchers do not directly
open a graphics device or save a plot. External pipeline diagnostics were not
treated as direct outputs of the archived launch scripts.

The 42 included files were confirmed by exact filename, output directory,
format, dimensions, file size, and checksum. No generated figure was assessed
for biological relevance, and no manuscript panel assignment is inferred.
