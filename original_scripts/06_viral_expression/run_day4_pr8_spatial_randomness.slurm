#!/bin/bash
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=31-00:00:00
#SBATCH --job-name=day4_pr8_spatial_test
#SBATCH --account=punim2394
#SBATCH --output=/data/projects/punim2394/lung_curio_mouse_pr8/viralgeneassess/%x_%j.log
#SBATCH --error=/data/projects/punim2394/lung_curio_mouse_pr8/viralgeneassess/%x_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this script with sbatch." >&2
  exit 1
fi

project_dir=/data/projects/punim2394/lung_curio_mouse_pr8
analysis_dir=${project_dir}/viralgeneassess
env_prefix=${project_dir}/OlliesFluRevisionCode/RCTD/conda_envs/rctd_spacexr_2.2
r_script=${analysis_dir}/test_day4_pr8_spatial_randomness.R

module purge
module load Anaconda3/2024.02-1

if [[ ! -x "${env_prefix}/bin/Rscript" ]]; then
  echo "Rscript is missing from the expected environment: ${env_prefix}" >&2
  exit 1
fi
if [[ ! -f "${r_script}" ]]; then
  echo "Analysis script is missing: ${r_script}" >&2
  exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd "${analysis_dir}"

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Project: ${project_dir}"
echo "R script: ${r_script}"
echo "Parameters: permutations=1999, k-neighbors=6, library-strata=10, minimum-group-size=100, seed=20260803"
echo "Inputs:"
printf '  %s\n' \
  "${project_dir}/processed_data/annotated_lungST_2_2_QCed.rds" \
  "${project_dir}/processed_data/annotated_lungST_2_3_QCed.rds"

"${env_prefix}/bin/Rscript" --vanilla -e '
cat(R.version.string, "\n")
cat("Seurat ", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("SeuratObject ", as.character(packageVersion("SeuratObject")), "\n", sep = "")
cat("Matrix ", as.character(packageVersion("Matrix")), "\n", sep = "")
cat("RANN ", as.character(packageVersion("RANN")), "\n", sep = "")
cat("ggplot2 ", as.character(packageVersion("ggplot2")), "\n", sep = "")
'

"${env_prefix}/bin/Rscript" --vanilla "${r_script}" \
  --n-permutations=1999 \
  --k-neighbors=6 \
  --library-strata=10 \
  --minimum-group-size=100 \
  --seed=20260803

echo "Finished: $(date --iso-8601=seconds)"
if command -v my-job-stats >/dev/null 2>&1; then
  my-job-stats -a -n -s || true
fi
