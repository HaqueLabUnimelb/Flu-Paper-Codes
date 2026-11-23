#!/bin/bash
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=31-00:00:00
#SBATCH --job-name=day10_rctd
#SBATCH --account=punim2394
#SBATCH --output=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD/RCTD_day10/logs/%x_%j.log
#SBATCH --error=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD/RCTD_day10/logs/%x_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this script with sbatch." >&2
  exit 1
fi

project_dir=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode
env_prefix=${project_dir}/RCTD/conda_envs/rctd_spacexr_2.2
r_script=${project_dir}/RCTD/RCTD_flu_d10_multi_spacexr2.2.R

module purge
module load Anaconda3/2024.02-1

if [[ ! -x "${env_prefix}/bin/Rscript" ]]; then
  echo "RCTD environment is missing: ${env_prefix}" >&2
  echo "Run RCTD/setup_spacexr2.2.sh first." >&2
  exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd "${project_dir}/RCTD"

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "R script: ${r_script}"
"${env_prefix}/bin/Rscript" -e '
stopifnot(getRversion() >= "4.3.0", getRversion() < "4.4.0")
stopifnot(packageVersion("spacexr") == "2.2.0")
stopifnot(packageVersion("Seurat") == "5.1.0")
cat(R.version.string, "\n")
cat("Seurat", as.character(packageVersion("Seurat")), "\n")
cat("spacexr", as.character(packageVersion("spacexr")), "\n")
'

"${env_prefix}/bin/Rscript" "${r_script}"

echo "Finished: $(date --iso-8601=seconds)"
if command -v my-job-stats >/dev/null 2>&1; then
  my-job-stats -a -n -s || true
fi
