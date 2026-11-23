# Running on a SLURM cluster

Submit from the root of the repository. The wrappers find the repository from the directory of submission; to
submit from elsewhere, export `FLU_REPO_DIR` first (`export FLU_REPO_DIR=/path/to/Flu-Paper-Codes`). SLURM
writes the logs of the R wrappers to `logs/` in the directory of submission, so create it there; the Curio
Seeker launchers write `slurm-<job>.out` there.

```bash
cp config/paths.example.env config/paths.env     # then edit the paths
mkdir -p logs

sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/01_preprocessing/lung_curioseeker_qc.sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/03_rctd/mouseDay4RCTD_full_v2.sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/03_rctd/mouseDay10RCTD_full_v3.sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/04_viral_mapping/lung_curioseeker_HKx31genomemapping.sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/04_viral_mapping/lung_curioseeker_SARSCOV2genomemapping.sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/05_revision_analyses/rctd/plot_individual_cell_types.sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/05_revision_analyses/rctd/run_same_input_doublet.sh
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION scripts/05_revision_analyses/rctd/compare_full_vs_doublet.sh
```

Resources. The wrappers request what the original runs requested (RCTD: one CPU, 300 GB, 31 days; the Curio
Seeker launchers: one CPU and 8 GB for preprocessing, eight CPUs and 50 GB for the viral references); override
these on the `sbatch` command line. Used by the original runs:

| Job | Run time | Peak memory |
|---|---|---|
| RCTD day 4 (two replicates, full mode) | 39 min | 4.5 GiB |
| RCTD day 10 | 15 min | 7.8 GiB |
| cell-type maps | 25 s | 0.5 GiB |
| Curio Seeker preprocessing, four samples (the launcher job; Nextflow submits the pipeline steps as jobs of their own) | 2 h 32 min | 1 GB |

The doublet-mode fit of the day-4 input takes about 4 hours on one CPU and about 3 GB.

Curio Seeker. Preprocessing used Curio Seeker v2.0.0 with `curioseeker_slurm_lung.config`, the viral references
v3.1.0 with `Bioinfo_Seeker_example_slurm.config` (both in `original_scripts/`); set `CURIOSEEKER_DIR` and
`CURIOSEEKER_NEXTFLOW_CONFIG` in `config/paths.env` for the step you run. The launchers need Nextflow and
Singularity or Apptainer on `PATH` (modules or a conda environment).

`RCTD_PREFLIGHT_ONLY=1` makes an RCTD script stop after reading its inputs and building the reference.

The R scripts can also be run directly, for example `Rscript --vanilla scripts/03_rctd/RCTD_flu_vis.R`, with
the variables of `config/paths.env` exported (`set -a; source config/paths.env; set +a`). Run with `Rscript`, the
feature-plot script writes the plots it only displays to `Rplots.pdf` (about 90 MB) in the working directory.
