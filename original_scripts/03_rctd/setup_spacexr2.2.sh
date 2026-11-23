#!/bin/bash
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=31-00:00:00
#SBATCH --job-name=rctd22_setup
#SBATCH --account=punim2394
#SBATCH --output=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD/logs/%x_%j.log
#SBATCH --error=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode/RCTD/logs/%x_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this script with sbatch." >&2
  exit 1
fi

project_dir=/data/projects/punim2394/lung_curio_mouse_pr8/OlliesFluRevisionCode
env_prefix=${project_dir}/RCTD/conda_envs/rctd_spacexr_2.2
package_cache=${project_dir}/RCTD/conda_pkgs
source_archive=${project_dir}/RCTD/software_source/spacexr-205c9acdcd5b9b794ae0521ab50473893e9e887f.tar.gz
environment_file=${project_dir}/RCTD/environment_spacexr2.2.yml
condarc_file=${project_dir}/RCTD/condarc_spacexr2.2.yml
expected_sha256=14da46d2dd14a469a4e548deff5319bdafbf23a4df3f825a01918af3262e355a

module purge
module load Anaconda3/2024.02-1

export CONDARC=${condarc_file}
export CONDA_PKGS_DIRS=${package_cache}
mkdir -p "${env_prefix%/*}" "${package_cache}"

echo "Started: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Conda: $(conda --version)"
echo "Environment: ${env_prefix}"

actual_sha256=$(sha256sum "${source_archive}" | awk '{print $1}')
if [[ "${actual_sha256}" != "${expected_sha256}" ]]; then
  echo "Unexpected spacexr source archive checksum: ${actual_sha256}" >&2
  exit 1
fi

if [[ -x "${env_prefix}/bin/R" ]]; then
  echo "Reusing existing conda environment: ${env_prefix}"
else
  conda env create \
    --prefix "${env_prefix}" \
    --file "${environment_file}" \
    --solver libmamba
fi

work_tmp=$(mktemp -d "/tmp/rctd22_${SLURM_JOB_ID}.XXXXXX")
source_dir=${work_tmp}/spacexr-205c9acdcd5b9b794ae0521ab50473893e9e887f
tar -xzf "${source_archive}" -C "${work_tmp}"
"${env_prefix}/bin/R" CMD INSTALL --no-multiarch "${source_dir}"

"${env_prefix}/bin/Rscript" -e '
stopifnot(getRversion() >= "4.3.0", getRversion() < "4.4.0")
stopifnot(packageVersion("spacexr") == "2.2.0")
stopifnot(packageVersion("Seurat") == "5.1.0")
cat(R.version.string, "\n")
cat("Seurat", as.character(packageVersion("Seurat")), "\n")
cat("spacexr", as.character(packageVersion("spacexr")), "\n")
'
echo "Running Seurat/spacexr input compatibility preflight"
RCTD_PREFLIGHT_ONLY=1 "${env_prefix}/bin/Rscript" \
  "${project_dir}/RCTD/RCTD_flu_d4_multi_spacexr2.2.R"
conda list --prefix "${env_prefix}" --explicit > "${project_dir}/RCTD/conda_spacexr2.2_explicit.txt"
echo "Finished: $(date --iso-8601=seconds)"
