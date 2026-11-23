# Differences between `scripts/` and `original_scripts/`

`original_scripts/` holds the files as they were run. `scripts/` holds the same analyses with the paths taken
from `config/paths.env` or from environment variables, so that they run outside the cluster they were written
on. No threshold, reference, sample grouping, RCTD setting, feature list, colour scale or figure size differs.
Every script of `scripts/` was compared line by line with its original; the R scripts and their wrappers were run
on the original data, the Curio Seeker launchers were checked by their command line
([`REPRODUCIBILITY.md`](REPRODUCIBILITY.md)).

| `scripts/` | Original | Differences |
|---|---|---|
| `01_preprocessing/lung_curioseeker_qc.sh` | `01_preprocessing/lung_curioseeker_qc.sh` | Paths from the configuration; account, partition and e-mail lines removed; `nextflow` must be on `PATH` (the original activates a conda environment); inputs are checked before the run. The `nextflow run` options are the same |
| `02_feature_plots/FeaturePlots_Flu.R` | `02_feature_plots/FeaturePlots_Flu_2026-07-01.R` | Input and output folders from `FLU_PROCESSED_DATA_DIR` and `FLU_FEATURE_PLOT_OUTPUT_DIR`. Six plots that are shown but not saved ask for a gene that is absent from the object; they are wrapped in `try()` so that `Rscript` continues as an interactive session does. The 42 saved plots are unchanged |
| `03_rctd/RCTD_flu_d4_multi_spacexr2.2.R`, `RCTD_flu_d10_multi_spacexr2.2.R` | same names in `03_rctd/` | Reference, input and output paths from environment variables; inputs are checked (files exist, two `SPATIAL` coordinates) |
| `03_rctd/mouseDay4RCTD_full_v2.sh`, `mouseDay10RCTD_full_v3.sh` | same names in `03_rctd/` | Account and absolute paths removed; the conda environment is given by `FLU_R_ENV_PREFIX`; the cluster's Anaconda module is loaded only when `FLU_ANACONDA_MODULE` is set; R is started with `--vanilla` |
| `03_rctd/RCTD_flu_vis.R` | `03_rctd/RCTD_flu_vis.R` | The original reads four results one after the other, so only the last one (day 10, full mode) is plotted; here one result is read (`FLU_RCTD_VIS_INPUT_RDS`, default the same one). The epithelial column is taken with either spelling (`Epithelial Cells` / `Epithelial cells`). The plot is saved when `FLU_RCTD_VIS_OUTPUT` is set |
| `04_viral_mapping/lung_curioseeker_HKx31genomemapping.sh`, `lung_curioseeker_SARSCOV2genomemapping.sh` | launchers in `04_viral_mapping/` | As for preprocessing: paths from the configuration, the same `nextflow run` options (the original SARS-CoV-2 run was started without `-resume`) and SLURM requests (eight CPUs, 50 GB, 31 days). One launcher per reference; the original runs used one launcher per run (X-31 day 4, X-31 day 10, SARS-CoV-2) and loaded Nextflow and Apptainer as modules |
| `05_revision_analyses/FishersExactTest.R` | `05_revision_analyses/FishersExactTest.R` | The three files are given as arguments or environment variables; the same calculation; the result is also written as a table, with the number of tested genes outside the background list |
| `05_revision_analyses/rctd/*.R` | same names in `05_revision_analyses/rctd/` | Paths from environment variables; the replicate labels and the objects are checked; default output folder of `plot_individual_cell_types.R` is `outputs/revision/rctd/day4_cell_type_plots` |
| `05_revision_analyses/rctd/*.sh`, `legacy/*.sh` | same names in `05_revision_analyses/rctd/` | As for the RCTD wrappers; they stop when not run under SLURM or with another R version than 4.3 |

In all scripts of `scripts/`:

- The R scripts find the repository from their own path and therefore run with `Rscript`, not with `source()`.
- The shell wrappers take the repository from `FLU_REPO_DIR`, from their own location when run with `bash`, or
  from the directory of submission when submitted with `sbatch`.
- A header says what the script does and how to run it; variables are quoted and `set -euo pipefail` is set.

The analyses of `original_scripts/06_viral_expression/` and the second version of the feature plots
(`original_scripts/02_feature_plots/FeaturePlots_Flu.R`) are archived as run only.
