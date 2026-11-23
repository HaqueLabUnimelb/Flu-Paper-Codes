# Software

## RCTD (steps 3, 5 and 6)

Conda environment of the original runs (built on 28 July 2026 by
`original_scripts/03_rctd/setup_spacexr2.2.sh`):

| Software | Version |
|---|---|
| R | 4.3.3 |
| Seurat | 5.1.0 |
| SeuratObject | 5.2.0 |
| spacexr | 2.2.0 |
| Matrix | 1.6-5 |
| ggplot2 | 3.5.2 |
| data.table | 1.17.8 |

- [`environment.yml`](environment.yml) is the environment file of that build;
  `original_scripts/03_rctd/conda_spacexr2.2_explicit.txt` lists every conda package with its build, and
  [`sessionInfo.txt`](sessionInfo.txt) is the session record of the environment.
- spacexr 2.2.0 is installed into the environment from source: commit
  `205c9acdcd5b9b794ae0521ab50473893e9e887f` of <https://github.com/dmcable/spacexr> ("v2.2.0"), with
  `R CMD INSTALL`. SHA-256 of the source archive used:
  `14da46d2dd14a469a4e548deff5319bdafbf23a4df3f825a01918af3262e355a`.

```bash
# the exact packages and builds of the original environment (linux-64)
conda create --prefix environment/rctd_spacexr_2.2 --file original_scripts/03_rctd/conda_spacexr2.2_explicit.txt
# or the same packages solved again: conda env create --prefix environment/rctd_spacexr_2.2 --file environment/environment.yml
curl -L -o spacexr.tar.gz https://github.com/dmcable/spacexr/archive/205c9acdcd5b9b794ae0521ab50473893e9e887f.tar.gz
tar -xzf spacexr.tar.gz
environment/rctd_spacexr_2.2/bin/R CMD INSTALL --no-multiarch spacexr-205c9acdcd5b9b794ae0521ab50473893e9e887f
```

The wrappers of `scripts/03_rctd/` and `scripts/05_revision_analyses/rctd/` take the environment from
`FLU_R_ENV_PREFIX` and check R 4.3, Seurat 5.1.0 (RCTD fits) and spacexr 2.2.0 before they run. The PR8
transcript scripts of `original_scripts/06_viral_expression/` ran in the same environment.

## Feature plots and RCTD overlay (step 2, `RCTD_flu_vis.R`)

These scripts load Seurat, spacexr, tidyverse, data.table and mltools. The 42 archived feature plots are
reproduced byte for byte with:

| Software | Version |
|---|---|
| R | 4.3.1 |
| Seurat | 4.3.0.1 |
| SeuratObject | 4.1.3 |
| ggplot2 | 3.5.0 |
| tidyverse | 2.0.0 |
| data.table | 1.14.8 |
| mltools | 0.3.5 |
| patchwork | 1.2.0 |
| spacexr | 2.2.1 |

With newer versions (R 4.5.0, Seurat 5.5.1, ggplot2 4.0.3) the script gives the same maps with larger points.

## Curio Seeker (steps 1 and 4)

| Step | Pipeline | Nextflow | Containers | STAR index |
|---|---|---|---|---|
| Preprocessing | Curio Seeker v2.0.0 | 23.10.1 | Singularity | 2.6.1d |
| Mapping to X-31 and SARS-CoV-2 references | Curio Seeker v3.1.0 | 23.04.2 | Apptainer 1.4.4 | 2.6.1d |

All versions with the script that uses them: [`DEPENDENCIES.tsv`](DEPENDENCIES.tsv).
