# Revision analyses

This module archives analyses added during manuscript revision. It is separated
from the primary preprocessing, feature-plot, RCTD-fitting and viral-mapping
workflows because these scripts answer reviewer-specific follow-up questions.

## Contents

| Path | Purpose |
|---|---|
| `FishersExactTest.R` | One-sided Fisher's exact test for overlap between two gene lists |
| `rctd/plot_individual_cell_types.R` | Spatial map for every cell-type proportion in each day-4 RCTD replicate |
| `rctd/run_same_input_doublet.R` | Rerun the fitted day-4 full-mode object in doublet mode without changing its inputs |
| `rctd/compare_full_vs_doublet.R` | Quantify agreement between the same-input full-mode and doublet-mode results |
| `rctd/legacy/` | Exploratory comparison of non-identically generated objects; retained with explicit limitations |

RDS inputs, generated tables, plots and logs are intentionally excluded from
Git. See [Input data requirements](../../docs/INPUT_DATA_REQUIREMENTS.md).

## Gene-list overlap test

The test asks whether the overlap between two unique gene sets is greater than
expected under random sampling from the supplied background universe. The first
column of each input text file is used. Empty identifiers and duplicates are
removed, and the script stops if either tested list contains a gene absent from
the background.

Run directly:

```bash
Rscript --vanilla scripts/05_revision_analyses/FishersExactTest.R \
  /path/to/background_genes.txt \
  /path/to/gene_list_1.txt \
  /path/to/gene_list_2.txt \
  outputs/revision/fisher_gene_overlap.tsv
```

Alternatively, set `FLU_FISHER_UNIVERSE`, `FLU_FISHER_GENE_LIST_1`,
`FLU_FISHER_GENE_LIST_2` and optionally `FLU_FISHER_OUTPUT` in
`config/paths.env`. The output records the universe size, both list sizes,
overlap size, odds ratio and one-sided p-value.

The background must represent the genes that could have entered either tested
list. It should not be replaced with all annotated genes unless that was the
actual selection universe.

## RCTD follow-up analyses

See [rctd/README.md](rctd/README.md) for execution order, SLURM resources and
interpretation limits. These scripts require the same R 4.3.x and spacexr 2.2.0
environment as the primary RCTD workflow.

No full analysis was rerun while preparing this public module.
