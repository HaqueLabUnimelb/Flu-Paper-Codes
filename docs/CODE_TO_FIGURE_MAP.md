# Code-to-figure map

No manuscript figure numbers were present in the supplied scripts, and no
figure-to-panel mapping was provided. The table therefore records only direct,
code-supported output relationships.

| Code section | Generated output | Repository status | Manuscript mapping |
|---|---|---|---|
| `FeaturePlots_Flu.R`: spatial QC | `mouse1_umis.jpeg`, `mouse2_umis.jpeg`, `mouse4_umis.jpeg`, `mouse5_umis.jpeg` | Expected; files not included | Unresolved |
| `FeaturePlots_Flu.R`: cytokine feature plots | Saved `mouse*_GENE.jpeg` files for the explicitly exported cytokine plots | Expected; files not included | Unresolved |
| `FeaturePlots_Flu.R`: B-cell-associated features | Saved `mouse*_Igkc.jpeg`, `mouse*_Ighm.jpeg`, and selected `mouse*_Ighg2c.jpeg` files | Expected; files not included | Unresolved |
| `RCTD_flu_vis.R`: cell-type overlay | Interactive plot; optional configured output | Not produced during repository preparation | Unresolved |

Curio Seeker, RCTD model-fitting, and viral-mapping scripts do not directly
open a graphics device or save a plot. Pipeline outputs may contain diagnostics,
but those are external workflow products and were not mapped here.

No generated figure was inspected for biological relevance, and no manuscript
panel assignment is inferred.
