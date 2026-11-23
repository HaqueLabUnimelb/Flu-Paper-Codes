# Input data

No data are in this repository. Put the inputs in `data/processed/` under the names below, or give their
paths in `config/paths.env`.

| File | Used by | Where it is |
|---|---|---|
| `lungST_2_2_QCed_edit.rds`, `lungST_2_3_QCed_edit.rds` (day 4), `lungST_4_2_QCed_edit.rds`, `lungST_4_3_QCed_edit.rds` (day 10) | RCTD | GEO [GSE341918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341918), supplementary files of samples GSM9920041 to GSM9920044 (named `GSM9920041_lungST_2_2_QCed_edit.rds` and so on) |
| `lungST_2_2_QCed.rds`, `lungST_2_3_QCed.rds`, `lungST_4_2_QCed.rds`, `lungST_4_3_QCed.rds` | feature plots, violin plots | lab's storage: the same preprocessing run as the GEO objects, another bead filter |
| `PR8_ref.rds` | RCTD | lab's storage: single-cell reference made of cells of GEO GSE202325 and GSE184384 |
| `annotated_lungST_*_QCed.rds` | PR8 transcript scripts (`original_scripts/06_viral_expression/`) | lab's storage |
| background list and SpatialDE gene lists | `FishersExactTest.R` | lab's storage |
| `Day4mouseRCTDoutmin30_editv1` | comparison with the earlier RCTD result | lab's storage |

The FASTQ files of the four lungs are in the SRA records of GSE341918. Details of every object:
[`../docs/INPUT_DATA_REQUIREMENTS.md`](../docs/INPUT_DATA_REQUIREMENTS.md) and
[`../docs/PROVENANCE.md`](../docs/PROVENANCE.md).
