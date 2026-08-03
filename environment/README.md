# Software environment

## Scope

`environment.yml` is the environment specification recovered from the RCTD
workflow. It pins R 4.3.x and Seurat 5.1.0 and lists the RCTD dependencies that
were recorded by the original setup process. It does not cover Curio Seeker,
Nextflow containers, or every package loaded by the plotting scripts.

spacexr was installed separately from a source archive corresponding to version
2.2.0. The recorded source directory identified commit
`205c9acdcd5b9b794ae0521ab50473893e9e887f`, and the archive SHA-256 was:

```text
14da46d2dd14a469a4e548deff5319bdafbf23a4df3f825a01918af3262e355a
```

The source archive is not distributed here, and its authoritative download
location was not recorded in the supplied scripts. Obtain spacexr 2.2.0 from an
authoritative source and verify the version before running RCTD.

## Confirmed RCTD environment

The existing environment was queried on 3 August 2026:

| Software/package | Version | Status |
|---|---:|---|
| R | 4.3.3 | Confirmed |
| Seurat | 5.1.0 | Confirmed |
| SeuratObject | 5.2.0 | Confirmed |
| spacexr | 2.2.0 | Confirmed |
| Matrix | 1.6-5 | Confirmed |
| ggplot2 | 3.5.2 | Confirmed |
| data.table | 1.17.8 | Confirmed |

The feature-plot and RCTD-visualization scripts also load `tidyverse` and
`mltools`. Those packages were not installed in the inspected RCTD
environment, so those scripts were parsed but not executed. Their exact
historical versions remain unknown.

## Creating the RCTD environment

```bash
conda env create \
  --prefix environment/rctd_spacexr_2.2 \
  --file environment/environment.yml
```

Install and verify spacexr 2.2.0 separately. The SLURM wrappers enforce R 4.3.x,
Seurat 5.1.0, and spacexr 2.2.0 before running. An HPC module named
`Anaconda3/2024.02-1` was recorded in the original wrappers; this is
site-specific and can be changed with `FLU_ANACONDA_MODULE`.

See `DEPENDENCIES.tsv` for evidence and status, and `sessionInfo.txt` for a
sanitized session record from the confirmed RCTD environment.
