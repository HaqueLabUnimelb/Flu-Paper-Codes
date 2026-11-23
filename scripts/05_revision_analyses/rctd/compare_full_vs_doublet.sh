#!/usr/bin/env bash
# Compare same-input full-mode and doublet-mode day-4 RCTD results.

#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=31-00:00:00
#SBATCH --job-name=rctd_mode_compare
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this script with sbatch." >&2
  exit 1
fi

# The repository: FLU_REPO_DIR when set; else the checkout that holds this file (run as "bash <file>"); else,
# under sbatch, which runs a copy of this file, the directory the job was submitted from.
script_rel=scripts/05_revision_analyses/rctd
script_name=compare_full_vs_doublet.sh
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." 2>/dev/null && pwd || true)
if [[ -n "${FLU_REPO_DIR:-}" ]]; then
  repo_dir=$(cd -- "${FLU_REPO_DIR}" && pwd)
  if [[ ! -f "${repo_dir}/${script_rel}/${script_name}" ]]; then
    echo "FLU_REPO_DIR is not this repository: ${repo_dir}" >&2
    exit 1
  fi
elif [[ -n "${here}" && -f "${here}/${script_rel}/${script_name}" ]]; then
  repo_dir=${here}
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/${script_rel}/${script_name}" ]]; then
  repo_dir=$(cd -- "${SLURM_SUBMIT_DIR}" && pwd)
else
  echo "Repository not found: submit from the root of the repository or set FLU_REPO_DIR." >&2
  exit 1
fi
script_dir="${repo_dir}/${script_rel}"
config_file=${FLU_CONFIG:-"${repo_dir}/config/paths.env"}

if [[ -f "${config_file}" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "${config_file}"
  set +a
fi

env_prefix=${FLU_R_ENV_PREFIX:-"${repo_dir}/environment/rctd_spacexr_2.2"}
r_script="${script_dir}/compare_full_vs_doublet.R"

# The original runs loaded the cluster's Anaconda module (FLU_ANACONDA_MODULE in config/paths.env); the R
# environment runs without it.
if [[ -n "${FLU_ANACONDA_MODULE:-}" ]] && type module >/dev/null 2>&1; then
  module purge
  module load "${FLU_ANACONDA_MODULE}"
fi
if [[ ! -x "${env_prefix}/bin/Rscript" ]]; then
  echo "R environment is missing: ${env_prefix}" >&2
  exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mkdir -p "${repo_dir}/outputs/revision/rctd"
cd "${repo_dir}"

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "R script: ${r_script}"
"${env_prefix}/bin/Rscript" --vanilla -e '
stopifnot(getRversion() >= "4.3.0", getRversion() < "4.4.0")
stopifnot(packageVersion("spacexr") == "2.2.0")
cat(R.version.string, "\n")
cat("spacexr", as.character(packageVersion("spacexr")), "\n")
'
"${env_prefix}/bin/Rscript" --vanilla "${r_script}"
echo "Finished: $(date --iso-8601=seconds)"
