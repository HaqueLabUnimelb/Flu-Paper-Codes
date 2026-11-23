#!/usr/bin/env bash
# Checks of this code archive that need no data: files, links, syntax, checksums. It changes nothing.
# Usage: [FLU_RSCRIPT=/path/to/Rscript] tests/validate_repository.sh
# Needs bash, find, grep, sha256sum, md5sum; git, R and Python 3 (with PyYAML) are used when present.
# Only the files of the repository are checked: in a git clone the tracked files and the files that are not
# ignored (so config/paths.env, data/, outputs/, logs/ and a conda environment under environment/ are left
# out); without git, every file except those folders.

set -uo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "${script_dir}/.." && pwd -P)
cd "${repo_dir}"

failures=0
printf 'check\tstatus\tevidence\n'

record() {
  printf '%s\t%s\t%s\n' "$1" "$2" "$3"
  [[ "$2" == "FAIL" ]] && failures=$((failures + 1))
  return 0
}

list_files() {
  if [[ "$(git -C "${repo_dir}" rev-parse --show-toplevel 2>/dev/null)" == "${repo_dir}" ]]; then
    git -C "${repo_dir}" ls-files -z --cached --others --exclude-standard
  else
    find . \( -path ./.git -o -path ./data -o -path ./outputs -o -path ./work -o -path ./logs \
      -o -path './environment/rctd_spacexr_2.2' -o -path ./config/paths.env -o -name Rplots.pdf \
      -o -name 'slurm-*.out' ! -path './original_scripts/*' \) -prune -o -type f -print0
    for f in data/README.md logs/validation_summary.txt; do
      [[ -f "${f}" ]] && printf '%s\0' "${f}"
    done
  fi
}
files=()
while IFS= read -r -d '' f; do
  f=${f#./}
  [[ -f "${f}" ]] && files+=("${f}")
done < <(list_files)

# pick <regex>: the files whose path matches the regex, NUL-separated
pick() {
  local f
  for f in "${files[@]}"; do
    [[ "${f}" =~ $1 ]] && printf '%s\0' "${f}"
  done
  return 0
}
binary_re='\.(jpe?g|png|pdf|gz|zip)$'
text_files() {
  local f
  for f in "${files[@]}"; do
    [[ "${f}" =~ ${binary_re} || "${f}" == "tests/validate_repository.sh" ]] || printf '%s\0' "${f}"
  done
  return 0
}

# 1. Required files
required_files=(
  README.md LICENSE CITATION.cff
  config/paths.example.env config/samples.example.tsv
  scripts/01_preprocessing/lung_curioseeker_qc.sh
  scripts/02_feature_plots/FeaturePlots_Flu.R
  scripts/03_rctd/RCTD_flu_d4_multi_spacexr2.2.R scripts/03_rctd/mouseDay4RCTD_full_v2.sh
  scripts/03_rctd/RCTD_flu_d10_multi_spacexr2.2.R scripts/03_rctd/mouseDay10RCTD_full_v3.sh
  scripts/03_rctd/RCTD_flu_vis.R
  scripts/04_viral_mapping/lung_curioseeker_HKx31genomemapping.sh
  scripts/04_viral_mapping/lung_curioseeker_SARSCOV2genomemapping.sh
  scripts/05_revision_analyses/FishersExactTest.R
  scripts/05_revision_analyses/rctd/plot_individual_cell_types.R
  scripts/05_revision_analyses/rctd/run_same_input_doublet.R
  scripts/05_revision_analyses/rctd/compare_full_vs_doublet.R
  scripts/05_revision_analyses/rctd/legacy/legacy_shared_spot_paired_ttest.R
  original_scripts/MANIFEST.tsv
  docs/PROVENANCE.md docs/CHANGES_FROM_ORIGINAL.md docs/FIGURE_MANIFEST.tsv
  figures/SHA256SUMS environment/environment.yml data/README.md
)
missing=()
for file in "${required_files[@]}"; do
  [[ -f "${file}" ]] || missing+=("${file}")
done
if (( ${#missing[@]} == 0 )); then
  record "Required files" "PASS" "${#required_files[@]} files exist"
else
  record "Required files" "FAIL" "missing: ${missing[*]}"
fi

# 2. Relative links in the Markdown files
broken=()
n_links=0
while IFS= read -r -d '' md; do
  while IFS= read -r link; do
    target=${link#*](}
    target=${target%%)*}
    target=${target%%#*}
    target=${target%% *}
    [[ -n "${target}" ]] || continue
    case "${target}" in http://*|https://*|mailto:*) continue ;; esac
    n_links=$((n_links + 1))
    [[ -e "$(dirname -- "${md}")/${target}" ]] || broken+=("${md} -> ${target}")
  done < <(grep -o '\]([^)]*)' "${md}" || true)
done < <(pick '\.md$')
if (( ${#broken[@]} == 0 )); then
  record "Markdown links" "PASS" "${n_links} relative links resolve"
else
  record "Markdown links" "FAIL" "${broken[*]}"
fi

# 3. Shell syntax (scripts/, tests/, config/ and the shell files of original_scripts/)
bad=()
n=0
while IFS= read -r -d '' f; do
  n=$((n + 1))
  bash -n "${f}" 2>/dev/null || bad+=("${f}")
done < <(pick '^(scripts|tests|original_scripts|config)/.*\.(sh|slurm|env)$')
if (( ${#bad[@]} == 0 )); then
  record "Shell syntax" "PASS" "${n} files pass bash -n"
else
  record "Shell syntax" "FAIL" "${bad[*]}"
fi

# 4. R syntax (scripts/ and original_scripts/)
rscript=${FLU_RSCRIPT:-$(command -v Rscript || true)}
if [[ -n "${rscript}" && -x "${rscript}" ]]; then
  bad=()
  n=0
  while IFS= read -r -d '' f; do
    n=$((n + 1))
    "${rscript}" --vanilla -e 'invisible(parse(file = commandArgs(TRUE)[1]))' "${f}" >/dev/null 2>&1 || bad+=("${f}")
  done < <(pick '^(scripts|original_scripts)/.*\.R$')
  if (( ${#bad[@]} == 0 )); then
    record "R syntax" "PASS" "${n} files parse with $("${rscript}" --version 2>&1 | head -n 1)"
  else
    record "R syntax" "FAIL" "${bad[*]}"
  fi
else
  record "R syntax" "SKIPPED" "no Rscript (set FLU_RSCRIPT)"
fi

# 5. Credentials and e-mail addresses
hits=$(text_files | xargs -0 -r grep -l -i -E \
  '(gh[pousr]_[A-Za-z0-9_]{20,}|github_pat_[A-Za-z0-9_]{20,}|glpat-[A-Za-z0-9_-]{20,}|xox[abposr]-[A-Za-z0-9-]{10,}|AKIA[0-9A-Z]{16}|BEGIN (RSA |OPENSSH |EC )?PRIVATE KEY|password[[:space:]]*=[[:space:]]*[^[:space:]]{6,})' \
  2>/dev/null || true)
if [[ -z "${hits}" ]]; then
  record "Credential patterns" "PASS" "none found"
else
  record "Credential patterns" "FAIL" "${hits//$'\n'/ }"
fi
hits=$(text_files | xargs -0 -r grep -l -E \
  '([A-Za-z0-9._%+-]+@([A-Za-z0-9-]+\.)+(com|org|net|edu|gov|io|ai|info|[a-z]{2})\.?([^A-Za-z0-9.-]|$)|--mail-user=[^<[:space:]])' \
  2>/dev/null | grep -v -x -e 'LICENSE' || true)
if [[ -z "${hits}" ]]; then
  record "E-mail addresses" "PASS" "none found"
else
  record "E-mail addresses" "FAIL" "${hits//$'\n'/ }"
fi

# 6. English text only
if command -v python3 >/dev/null 2>&1; then
  hits=$(text_files | python3 -c '
import re, sys
cjk = re.compile("[\u3040-\u30ff\u3400-\u9fff\uac00-\ud7af]")
for path in sys.stdin.buffer.read().decode("utf-8", "surrogateescape").split("\0"):
    if path:
        with open(path, encoding="utf-8", errors="ignore") as fh:
            if cjk.search(fh.read()):
                print(path)
')
  if [[ -z "${hits}" ]]; then
    record "English text" "PASS" "no CJK character in the text files"
  else
    record "English text" "FAIL" "${hits//$'\n'/ }"
  fi
else
  record "English text" "SKIPPED" "no Python 3"
fi

# 7. No research data, no oversized file
hits=$(pick '\.([Rr][Dd][Ss]|[Rr][Dd]ata|[Rr][Dd][Aa]|fastq|fastq\.gz|fq|fq\.gz|bam|h5|h5ad)$' | tr '\0' ' ')
if [[ -z "${hits}" ]]; then
  record "Data files" "PASS" "no data object in the repository"
else
  record "Data files" "FAIL" "${hits}"
fi
hits=""
for f in "${files[@]}"; do
  (( $(stat -c %s "${f}") > 99 * 1024 * 1024 )) && hits+="${f} "
done
if [[ -z "${hits}" ]]; then
  record "File size" "PASS" "no file above 99 MiB"
else
  record "File size" "FAIL" "${hits}"
fi

# 8. The configurable scripts carry no site path (original_scripts/ keeps the paths of the cluster it ran on)
hits=$(pick '^(scripts|config)/' | xargs -0 -r grep -l -i -E '(/data/(projects|gpfs|scratch)/|punim[0-9]{4})' \
  2>/dev/null || true)
if [[ -z "${hits}" ]]; then
  record "Site paths in scripts/ and config/" "PASS" "none found"
else
  record "Site paths in scripts/ and config/" "FAIL" "${hits//$'\n'/ }"
fi

# 9. Executable permissions
hits=""
while IFS= read -r -d '' f; do
  [[ -x "${f}" ]] || hits+="${f} "
done < <(pick '^(scripts|tests)/.*\.(sh|R)$')
if [[ -z "${hits}" ]]; then
  record "Executable permissions" "PASS" "all files of scripts/ and tests/"
else
  record "Executable permissions" "FAIL" "${hits}"
fi

# 10. Figures: manifest and checksums
missing=()
n=0
while IFS=$'\t' read -r file _rest; do
  [[ "${file}" == "file" ]] && continue
  n=$((n + 1))
  [[ -f "${file}" ]] || missing+=("${file}")
done < docs/FIGURE_MANIFEST.tsv
n_disk=$(pick '^figures/generated/' | tr '\0' '\n' | grep -c -v '/\.gitkeep$')
if (( ${#missing[@]} == 0 && n == n_disk )); then
  record "Figure manifest" "PASS" "${n} files listed, ${n_disk} in figures/generated"
else
  record "Figure manifest" "FAIL" "listed ${n}, in figures/generated ${n_disk}, missing: ${missing[*]:-none}"
fi
manifest_sums=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "sha256") c = i; next } { print $c "  " $1 }' \
  docs/FIGURE_MANIFEST.tsv | sort)
if ! sha256sum -c figures/SHA256SUMS >/dev/null 2>&1; then
  record "Figure checksums" "FAIL" "sha256sum -c figures/SHA256SUMS reports a difference"
elif [[ "${manifest_sums}" != "$(sort figures/SHA256SUMS)" ]]; then
  record "Figure checksums" "FAIL" "docs/FIGURE_MANIFEST.tsv and figures/SHA256SUMS list different files or checksums"
else
  record "Figure checksums" "PASS" "$(wc -l < figures/SHA256SUMS) files match figures/SHA256SUMS and the manifest"
fi

# 11. Original scripts: every file is listed and equals its recorded md5
bad=()
n=0
while IFS=$'\t' read -r path _source _owner _modified _size _md5_source md5_copy _note; do
  [[ "${path}" == "path" ]] && continue
  n=$((n + 1))
  if [[ ! -f "${path}" ]] || [[ "$(md5sum < "${path}" | cut -d ' ' -f 1)" != "${md5_copy}" ]]; then
    bad+=("${path}")
  fi
done < original_scripts/MANIFEST.tsv
n_disk=$(pick '^original_scripts/' | tr '\0' '\n' | grep -c -v -x -e 'original_scripts/MANIFEST.tsv' -e 'original_scripts/README.md')
if (( ${#bad[@]} == 0 && n == n_disk )); then
  record "Original scripts" "PASS" "${n} files equal the md5 in original_scripts/MANIFEST.tsv"
else
  record "Original scripts" "FAIL" "listed ${n}, in original_scripts ${n_disk}, different: ${bad[*]:-none}"
fi

# 12. YAML syntax
if command -v python3 >/dev/null 2>&1 && python3 -c 'import yaml' >/dev/null 2>&1; then
  if python3 -c 'import yaml; yaml.safe_load(open("environment/environment.yml")); yaml.safe_load(open("CITATION.cff"))' 2>/dev/null; then
    record "YAML syntax" "PASS" "environment.yml and CITATION.cff parse"
  else
    record "YAML syntax" "FAIL" "a YAML file does not parse"
  fi
else
  record "YAML syntax" "SKIPPED" "no Python with PyYAML"
fi

if (( failures > 0 )); then
  record "Overall" "FAIL" "${failures} failed check(s)"
  exit 1
fi
record "Overall" "PASS" "no failed check"
