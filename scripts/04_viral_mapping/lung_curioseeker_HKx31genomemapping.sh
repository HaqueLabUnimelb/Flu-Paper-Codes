#!/usr/bin/env bash
# Curio Seeker mapping against the HKx31 reference.
# Required variables are documented in config/paths.example.env.
# Submit with sbatch; this Nextflow workflow must not be run on a login node.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=31-00:00:00
#SBATCH --job-name=curioseeker_hkx31

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this script with sbatch." >&2
  exit 1
fi

# The repository: FLU_REPO_DIR when set; else the checkout that holds this file (run as "bash <file>"); else,
# under sbatch, which runs a copy of this file, the directory the job was submitted from.
script_rel=scripts/04_viral_mapping
script_name=lung_curioseeker_HKx31genomemapping.sh
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd || true)
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

require_var() {
  local name=$1
  if [[ -z "${!name:-}" ]]; then
    echo "Required variable is not set: ${name}. See config/paths.example.env." >&2
    exit 1
  fi
}

for name in CURIOSEEKER_DIR HKX31_SAMPLESHEET HKX31_OUTPUT_DIR \
  VIRAL_MAPPING_WORK_DIR VIRAL_IGENOMES_BASE CURIOSEEKER_NEXTFLOW_CONFIG; do
  require_var "${name}"
done

for required_file in "${CURIOSEEKER_DIR}/main.nf" "${HKX31_SAMPLESHEET}" \
  "${VIRAL_IGENOMES_BASE}" "${CURIOSEEKER_NEXTFLOW_CONFIG}"; do
  if [[ ! -e "${required_file}" ]]; then
    echo "Required input was not found: ${required_file}" >&2
    exit 1
  fi
done

if ! command -v nextflow >/dev/null 2>&1; then
  echo "nextflow is not available in PATH." >&2
  exit 1
fi

mkdir -p "${HKX31_OUTPUT_DIR}" "${VIRAL_MAPPING_WORK_DIR}"
cd "${CURIOSEEKER_DIR}"

nextflow run main.nf \
  --input "${HKX31_SAMPLESHEET}" \
  --outdir "${HKX31_OUTPUT_DIR}" \
  -work-dir "${VIRAL_MAPPING_WORK_DIR}" \
  --igenomes_base "${VIRAL_IGENOMES_BASE}" \
  -resume \
  -profile singularity \
  -config "${CURIOSEEKER_NEXTFLOW_CONFIG}"

if command -v my-job-stats >/dev/null 2>&1; then
  my-job-stats -a -n -s || true
fi
