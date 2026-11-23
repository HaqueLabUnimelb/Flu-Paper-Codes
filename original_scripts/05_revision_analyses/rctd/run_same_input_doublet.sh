#!/bin/bash
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=31-00:00:00
#SBATCH --job-name=rctd_same_doublet
#SBATCH --account=punim2394
#SBATCH --output=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD/logs/%x_%j.log
#SBATCH --error=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD/logs/%x_%j.log

set -euo pipefail

project_dir=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode
env_prefix=${project_dir}/RCTD/conda_envs/rctd_spacexr_2.2
r_script=${project_dir}/RCTD/run_same_input_doublet.R

module purge
module load Anaconda3/2024.02-1

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
"${env_prefix}/bin/Rscript" -e 'stopifnot(packageVersion("spacexr") == "2.2.0")'
"${env_prefix}/bin/Rscript" "${r_script}"
echo "Finished: $(date --iso-8601=seconds)"
