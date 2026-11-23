# Input data

No data are in this repository. Where each input is: [`../data/README.md`](../data/README.md); how it was made:
[`PROVENANCE.md`](PROVENANCE.md).

## Spatial Seurat objects

| Sample | Day | RCTD input (GEO GSE341918) | Feature-plot input |
|---|---|---|---|
| mouse1 | 4 | `lungST_2_2_QCed_edit.rds` | `lungST_2_2_QCed.rds` |
| mouse2 | 4 | `lungST_2_3_QCed_edit.rds` | `lungST_2_3_QCed.rds` |
| mouse4 | 10 | `lungST_4_2_QCed_edit.rds` | `lungST_4_2_QCed.rds` |
| mouse5 | 10 | `lungST_4_3_QCed_edit.rds` | `lungST_4_3_QCed.rds` |

Each object has an `RNA` assay with counts and a two-dimensional `SPATIAL` reduction. The two sets come from
the same preprocessing run and differ in the beads kept (bead numbers in [`PROVENANCE.md`](PROVENANCE.md)).
GEO adds the sample accession to the file names (`GSM9920041_lungST_2_2_QCed_edit.rds`); remove it, or give the
paths in `config/paths.env`.

## Single-cell reference

`PR8_ref.rds` is a Seurat object with RNA counts and the metadata columns `timepoint_ref` (`4dpi`, `10dpi`) and
`Annotation_relabel` (cell type). It holds 17,826 cells of GEO GSE202325 (3 days after infection) and
GSE184384 (5, 9 and 11 days after infection); composition in [`PROVENANCE.md`](PROVENANCE.md).

## Revision analyses

- Fisher's exact test: three text files with one gene per line (background list, list 1, list 2).
- RCTD cell-type maps and mode comparison: the day-4 result `Day4mouseRCTDoutmin30_multi.rds` of
  `scripts/03_rctd/`; the comparison also reads the doublet-mode result written by `run_same_input_doublet.R`.
- Comparison with the earlier result: the doublet-mode object `Day4mouseRCTDoutmin30_editv1`.

## Curio Seeker

Preprocessing and viral mapping need the Curio Seeker workflow (`main.nf`), a samplesheet, a Nextflow
configuration, the reference folder given as `--igenomes_base` (it holds one sub-folder per genome key of the
samplesheet) and the container image of the workflow. The samplesheets and configurations of the original runs
are in `original_scripts/01_preprocessing/` and `original_scripts/04_viral_mapping/`; the FASTQ files are in
the SRA records of GSE341918.
