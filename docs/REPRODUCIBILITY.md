# Reproducibility status

## Verified

- All nine supplied source scripts existed, were readable, and were copied
  without modifying the originals.
- Five Shell scripts passed `bash -n`.
- Four R scripts were parsed successfully with R 4.3.3.
- The inspected RCTD environment contained Seurat 5.1.0 and spacexr 2.2.0.
- Public scripts use configuration variables instead of institutional paths.
- Scientific parameters visible in the supplied scripts were retained.

## Inferred

- Curio Seeker version 3.0.0 is inferred from the recorded workflow directory.
- A Singularity-compatible runtime is inferred from the Nextflow profile.
- Sample/time-point labels are taken from comments and filenames in the code.

## Not tested

- No complete Nextflow, feature-plot, RCTD, or viral-mapping analysis was run.
- No generated figure was compared with a manuscript panel.
- `shellcheck` was unavailable in the current environment.
- The plotting scripts were not executed because their full package environment
  and input data were not available together.
- GitHub Actions or another remote CI service is not configured.

## Missing or unresolved

- Raw and processed data are not distributed here.
- Some inputs recorded by the source scripts were not found during the audit.
- Exact historical versions of all Curio Seeker and plotting dependencies are
  unavailable.
- The authoritative spacexr source-archive location is unresolved.
- Manuscript title, citation, data accessions, and figure panel mapping require
  author confirmation.

Exact reproduction may depend on unavailable inputs, reference resources,
containers, and historical software versions.
