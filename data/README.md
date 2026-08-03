# Input data

No research data are committed to this repository.

Place local inputs under `data/processed/` using the filenames below, or set
the corresponding variables in `config/paths.env`:

- `lungST_2_2_QCed.rds` and `lungST_2_3_QCed.rds` for day-4 feature plots.
- `lungST_4_2_QCed.rds` and `lungST_4_3_QCed.rds` for day-10 feature plots.
- Edited versions of those four files, with `_edit` before `.rds`, for RCTD.
- `PR8_ref.rds`, containing the single-cell reference and the metadata fields
  expected by the RCTD scripts.
- Curio Seeker samplesheets, workflow configuration, and compatible reference
  resources for preprocessing and viral mapping.

The scripts also require a `SPATIAL` reduction in each spatial Seurat object.
The RCTD reference is expected to contain `timepoint_ref` and
`Annotation_relabel` metadata. These requirements are inferred directly from
the code; the repository does not validate or redistribute the objects.

The associated dataset has been deposited in NCBI GEO under Series accession
[GSE341918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341918).
Users should confirm which required repository inputs are represented in that
deposit and whether the record has been publicly released. Raw FASTQ files,
alignment files, reference genomes, Seurat/RCTD objects, pipeline work
directories, and credentials are intentionally excluded from GitHub.
