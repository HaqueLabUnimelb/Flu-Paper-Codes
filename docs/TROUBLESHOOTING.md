# Troubleshooting

| Message or symptom | What to do |
|---|---|
| "Required variable is not set" | Copy `config/paths.example.env` to `config/paths.env`, fill in the paths and submit from the root of the repository; or set `FLU_CONFIG` to the configuration file |
| "Repository not found" | Submit from the root of the repository, or export `FLU_REPO_DIR` before `sbatch` |
| "Required RDS input(s) not found" | Check the file names against `config/samples.example.tsv`; the RCTD scripts read the `_QCed_edit.rds` objects, the feature plots the `_QCed.rds` objects |
| The version check of an RCTD wrapper stops the job | The wrappers require R 4.3, Seurat 5.1.0 and spacexr 2.2.0; point `FLU_R_ENV_PREFIX` to that environment |
| `module load` fails | `FLU_ANACONDA_MODULE` names a module of the Spartan cluster; remove it from `config/paths.env` on other systems |
| `library(tidyverse)` or `library(mltools)` fails | The feature-plot and visualisation scripts need them; the RCTD environment does not hold them (see `environment/README.md`) |
| "nextflow is not available in PATH" | Load Nextflow and Singularity or Apptainer before `sbatch` |
| SLURM cannot open the log file | Run `mkdir -p logs` in the directory you submit from (the root of the repository) |
