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

Public accession identifiers and controlled-access conditions were not
available during repository preparation. The manuscript authors must add the
confirmed data-availability statement. Raw FASTQ files, alignment files,
reference genomes, Seurat/RCTD objects, pipeline work directories, and
credentials are intentionally excluded from GitHub.
