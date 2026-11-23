# Comparison with the earlier doublet-mode result

`legacy_shared_spot_paired_ttest.R` compares the day-4 full-mode result with the earlier doublet-mode result
`Day4mouseRCTDoutmin30_editv1` (doublet mode, `UMI_min = 75`, 16,292 and 14,421 beads) on the beads that its fit
does not reject:
proportions of every cell type in both results and a paired t-test per cell type. Because input and settings
differ, a difference cannot be assigned to the mode alone; [`../compare_full_vs_doublet.R`](../compare_full_vs_doublet.R)
compares the modes on the same input. Original run: SLURM job 28266050, 29 July 2026
([`../../../../original_scripts/05_revision_analyses/rctd/`](../../../../original_scripts/05_revision_analyses/rctd)).
