# Input data requirements

## Spatial Seurat objects

The feature-plot script expects:

| Sample label | Time point | Filename |
|---|---|---|
| mouse1 | day 4 | `lungST_2_2_QCed.rds` |
| mouse2 | day 4 | `lungST_2_3_QCed.rds` |
| mouse4 | day 10 | `lungST_4_2_QCed.rds` |
| mouse5 | day 10 | `lungST_4_3_QCed.rds` |

The RCTD scripts expect corresponding files with `_edit` before the extension.
Each object must contain an RNA assay and a two-dimensional `SPATIAL`
reduction. These requirements are derived from code access, not from a
biological review of the objects.

## Single-cell reference

`PR8_ref.rds` must be a Seurat object containing:

- `timepoint_ref`, including `4dpi` and `10dpi`;
- `Annotation_relabel`, used as the cell-type annotation;
- RNA count data compatible with the spatial objects.

The reference object is not included.

## Revision-analysis inputs

The Fisher overlap script expects three text files:

- a background universe containing every gene that could have been selected;
- gene list 1;
- gene list 2.

The first column is read without a header. Identifiers are trimmed and
deduplicated. Every gene in either tested list must occur in the universe.

The primary RCTD revision workflow expects the day-4 full-mode
`RCTD.replicates` object produced by `scripts/03_rctd/`. The doublet rerun
creates a second object from those same fitted inputs. The comparison script
requires both objects and assumes identical replicate and spot ordering.

The optional legacy comparison additionally requires an older day-4
doublet-mode object. That object was generated from a non-identical spot set or
preprocessing path and is not included. The legacy workflow is provided for
provenance, not as the preferred mode-sensitivity analysis.

## Curio Seeker inputs

Preprocessing and viral mapping require:

- an external Curio Seeker workflow containing `main.nf`;
- a workflow-compatible CSV samplesheet;
- a Nextflow configuration file;
- the relevant host or viral reference resources;
- container images retrievable by the configured Singularity profile.

The meaning expected by Curio Seeker's `--igenomes_base` option could not be
verified from the supplied launch scripts alone. Configure it according to the
matching Curio Seeker version.

## Availability

Four QCed feature-plot RDS inputs and the single-cell reference existed during
the technical audit. The four edited RCTD spatial inputs were not found at the
paths recorded in the supplied scripts, although existing day-4 and day-10 RCTD
result objects were found. The viral-mapping samplesheets recorded by the
scripts were not found.

The associated dataset has been deposited in NCBI GEO under Series accession
[GSE341918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341918).
Users should confirm which required local inputs are represented in the deposit
and whether the record has been publicly released. No data files are committed
to GitHub.
