# Revision analyses

Analyses added during the revision of the manuscript. Originals as run:
[`../../original_scripts/05_revision_analyses/`](../../original_scripts/05_revision_analyses).

| Path | Purpose |
|---|---|
| `FishersExactTest.R` | One-sided Fisher's exact test of the overlap of two gene lists |
| `rctd/plot_individual_cell_types.R` | Map of every cell-type proportion in each day-4 RCTD replicate |
| `rctd/run_same_input_doublet.R` | Fits the day-4 full-mode object again in doublet mode, on the same input |
| `rctd/compare_full_vs_doublet.R` | Agreement between the full-mode and the doublet-mode result of the same input |
| `rctd/legacy/` | Comparison of the full-mode result with an earlier doublet-mode result |

## Gene-list overlap test

```bash
Rscript --vanilla scripts/05_revision_analyses/FishersExactTest.R \
  /path/to/background_genes.txt /path/to/gene_list_1.txt /path/to/gene_list_2.txt \
  outputs/revision/fisher_gene_overlap.tsv
```

or set `FLU_FISHER_UNIVERSE`, `FLU_FISHER_GENE_LIST_1`, `FLU_FISHER_GENE_LIST_2` and `FLU_FISHER_OUTPUT`.
The first column of each file is read. With N genes in the background list, a and b genes in the two lists and
k genes in both, the test is `fisher.test(matrix(c(k, a - k, b - k, N - (a + b - k)), nrow = 2),
alternative = "greater")`. The script prints the test and writes N, a, b, k, the number of tested genes that
are not in the background list, the odds ratio and the p-value.

On the original inputs (background list of 21,214 genes, SpatialDE gene lists of 61 and 48 genes, six genes of
each not in the background list) the script gives an overlap of 29 genes, odds ratio 1,024, p = 7.4e-65.

## RCTD analyses

See [`rctd/README.md`](rctd/README.md). They need R 4.3 and spacexr 2.2.0, as the RCTD fits
([`../../environment/README.md`](../../environment/README.md)).
