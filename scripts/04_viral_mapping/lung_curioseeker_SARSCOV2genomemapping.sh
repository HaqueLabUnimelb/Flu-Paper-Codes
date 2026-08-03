#!/usr/bin/env bash
# Curio Seeker mapping against the SARS-CoV-2 reference.
# Required variables are documented in config/paths.example.env.
# Submit with sbatch; this Nextflow workflow must not be run on a login node.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192
#SBATCH --time=14-00:00:00
#SBATCH --job-name=curioseeker_sarscov2

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

require_var() {
  local name=$1
  if [[ -z "${!name:-}" ]]; then
    echo "Required variable is not set: ${name}. See config/paths.example.env." >&2
    exit 1
  fi
}

for name in CURIOSEEKER_DIR SARSCOV2_SAMPLESHEET SARSCOV2_OUTPUT_DIR \
  VIRAL_MAPPING_WORK_DIR VIRAL_IGENOMES_BASE CURIOSEEKER_NEXTFLOW_CONFIG; do
  require_var "${name}"
done

for required_file in "${CURIOSEEKER_DIR}/main.nf" "${SARSCOV2_SAMPLESHEET}" \
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

mkdir -p "${SARSCOV2_OUTPUT_DIR}" "${VIRAL_MAPPING_WORK_DIR}"
cd "${CURIOSEEKER_DIR}"

nextflow run main.nf \
  --input "${SARSCOV2_SAMPLESHEET}" \
  --outdir "${SARSCOV2_OUTPUT_DIR}" \
  -work-dir "${VIRAL_MAPPING_WORK_DIR}" \
  --igenomes_base "${VIRAL_IGENOMES_BASE}" \
  -resume \
  -profile singularity \
  -config "${CURIOSEEKER_NEXTFLOW_CONFIG}"

if command -v my-job-stats >/dev/null 2>&1; then
  my-job-stats -a -n -s || true
fi
