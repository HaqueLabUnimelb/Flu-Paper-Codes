#!/usr/bin/env bash
# Run day-4 RCTD in full doublet mode. This is a 300 GB SLURM job.

#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=31-00:00:00
#SBATCH --job-name=day4_rctd22
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this script with sbatch." >&2
  exit 1
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "${script_dir}/../.." && pwd)
config_file=${FLU_CONFIG:-"${repo_dir}/config/paths.env"}

if [[ -f "${config_file}" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "${config_file}"
  set +a
fi

env_prefix=${FLU_R_ENV_PREFIX:-"${repo_dir}/environment/rctd_spacexr_2.2"}
r_script="${script_dir}/RCTD_flu_d4_multi_spacexr2.2.R"

if type module >/dev/null 2>&1; then
  module purge
  module load "${FLU_ANACONDA_MODULE:-Anaconda3/2024.02-1}"
fi

if [[ ! -x "${env_prefix}/bin/Rscript" ]]; then
  echo "R environment is missing: ${env_prefix}" >&2
  echo "Set FLU_R_ENV_PREFIX to an environment containing R 4.3, Seurat 5.1.0 and spacexr 2.2.0." >&2
  exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mkdir -p "${repo_dir}/logs" "${FLU_RCTD_OUTPUT_DIR:-${repo_dir}/outputs/rctd}"
cd "${repo_dir}"

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "R script: ${r_script}"
"${env_prefix}/bin/Rscript" --vanilla -e '
stopifnot(getRversion() >= "4.3.0", getRversion() < "4.4.0")
stopifnot(packageVersion("spacexr") == "2.2.0")
stopifnot(packageVersion("Seurat") == "5.1.0")
cat(R.version.string, "\n")
cat("Seurat", as.character(packageVersion("Seurat")), "\n")
cat("spacexr", as.character(packageVersion("spacexr")), "\n")
'
"${env_prefix}/bin/Rscript" --vanilla "${r_script}"

echo "Finished: $(date --iso-8601=seconds)"
if command -v my-job-stats >/dev/null 2>&1; then
  my-job-stats -a -n -s || true
fi
