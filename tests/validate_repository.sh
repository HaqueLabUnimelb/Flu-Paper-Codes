#!/usr/bin/env bash
# Lightweight, non-destructive validation for this public code archive.

set -uo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "${script_dir}/.." && pwd)
cd "${repo_dir}"

failures=0
warnings=0

printf 'check\tstatus\tevidence\taction_taken\tremaining_issue\n'

record() {
  local check=$1
  local status=$2
  local evidence=$3
  local action=$4
  local issue=$5
  printf "%s\t%s\t%s\t%s\t%s\n" "${check}" "${status}" "${evidence}" "${action}" "${issue}"
  case "${status}" in
    FAIL) failures=$((failures + 1)) ;;
    WARNING) warnings=$((warnings + 1)) ;;
  esac
}

required_files=(
  README.md
  LICENSE
  CITATION.cff
  PROMPTS/CODEX_REPRODUCTION_PROMPT.md
  config/paths.example.env
  config/samples.example.tsv
  scripts/01_preprocessing/lung_curioseeker_qc.sh
  scripts/02_feature_plots/FeaturePlots_Flu.R
  scripts/03_rctd/RCTD_flu_d4_multi_spacexr2.2.R
  scripts/03_rctd/mouseDay4RCTD_full_v2.sh
  scripts/03_rctd/RCTD_flu_d10_multi_spacexr2.2.R
  scripts/03_rctd/mouseDay10RCTD_full_v3.sh
  scripts/03_rctd/RCTD_flu_vis.R
  scripts/04_viral_mapping/lung_curioseeker_HKx31genomemapping.sh
  scripts/04_viral_mapping/lung_curioseeker_SARSCOV2genomemapping.sh
  docs/CODE_TO_FIGURE_MAP.md
  docs/FIGURE_MANIFEST.tsv
  figures/SHA256SUMS
  environment/environment.yml
  data/README.md
)

missing_required=()
for file in "${required_files[@]}"; do
  [[ -f "${file}" ]] || missing_required+=("${file}")
done
if (( ${#missing_required[@]} == 0 )); then
  record "Required repository files" "PASS" "All ${#required_files[@]} required files exist" "None" "None"
else
  record "Required repository files" "FAIL" "${#missing_required[@]} files are missing" "None" "${missing_required[*]}"
fi

readme_link=$(sed -n '3p' README.md)
if [[ "${readme_link}" == '**AI-assisted reproduction:** [Open the Codex reproduction prompt](PROMPTS/CODEX_REPRODUCTION_PROMPT.md)' ]]; then
  record "README prompt placement" "PASS" "Prompt link is the first content after the title" "None" "None"
else
  record "README prompt placement" "FAIL" "README line 3 does not match the required link" "None" "Move the prompt link directly below the title"
fi

broken_links=()
while IFS= read -r markdown_file; do
  while IFS= read -r raw_link; do
    target=${raw_link#](}
    target=${target%%#*}
    target=${target%% *}
    [[ -n "${target}" ]] || continue
    case "${target}" in
      http://*|https://*|mailto:*|\#*) continue ;;
    esac
    link_base=$(dirname -- "${markdown_file}")
    if [[ ! -e "${link_base}/${target}" ]]; then
      broken_links+=("${markdown_file} -> ${target}")
    fi
  done < <(rg --no-filename --only-matching '\]\([^)]+' "${markdown_file}" || true)
done < <(find . -path './.git' -prune -o -name '*.md' -type f -print)
if (( ${#broken_links[@]} == 0 )); then
  record "Internal Markdown links" "PASS" "All relative Markdown links resolve" "None" "None"
else
  record "Internal Markdown links" "FAIL" "${#broken_links[@]} broken link(s)" "None" "${broken_links[*]}"
fi

if bash -n config/paths.example.env; then
  record "Configuration syntax" "PASS" "config/paths.example.env passes bash -n" "None" "None"
else
  record "Configuration syntax" "FAIL" "Example environment configuration is invalid" "None" "Correct shell syntax"
fi

shell_failures=()
while IFS= read -r -d '' shell_file; do
  bash -n "${shell_file}" || shell_failures+=("${shell_file}")
done < <(find scripts tests -name '*.sh' -type f -print0)
if (( ${#shell_failures[@]} == 0 )); then
  record "Shell syntax" "PASS" "All public Shell scripts pass bash -n" "None" "None"
else
  record "Shell syntax" "FAIL" "${#shell_failures[@]} Shell file(s) failed" "None" "${shell_failures[*]}"
fi

rscript_path=${FLU_RSCRIPT:-}
if [[ -z "${rscript_path}" ]] && command -v Rscript >/dev/null 2>&1; then
  rscript_path=$(command -v Rscript)
fi
if [[ -n "${rscript_path}" && -x "${rscript_path}" ]]; then
  r_failures=()
  while IFS= read -r -d '' r_file; do
    "${rscript_path}" --vanilla -e 'parse(file=commandArgs(TRUE)[1])' "${r_file}" >/dev/null ||
      r_failures+=("${r_file}")
  done < <(find scripts -name '*.R' -type f -print0)
  if (( ${#r_failures[@]} == 0 )); then
    record "R syntax" "PASS" "All public R scripts parsed with ${rscript_path}" "None" "None"
  else
    record "R syntax" "FAIL" "${#r_failures[@]} R file(s) failed to parse" "None" "${r_failures[*]}"
  fi
else
  record "R syntax" "NOT TESTED" "No executable Rscript was supplied or found" "None" "Set FLU_RSCRIPT and rerun"
fi

if command -v shellcheck >/dev/null 2>&1; then
  shellcheck_failures=()
  while IFS= read -r -d '' shell_file; do
    shellcheck "${shell_file}" || shellcheck_failures+=("${shell_file}")
  done < <(find scripts tests -name '*.sh' -type f -print0)
  if (( ${#shellcheck_failures[@]} == 0 )); then
    record "ShellCheck" "PASS" "All Shell scripts pass shellcheck" "None" "None"
  else
    record "ShellCheck" "WARNING" "${#shellcheck_failures[@]} Shell file(s) reported findings" "None" "${shellcheck_failures[*]}"
  fi
else
  record "ShellCheck" "NOT TESTED" "shellcheck is not installed" "None" "Install shellcheck in an authorized environment to test"
fi

secret_matches=$(rg -l -i "(gh[pousr]_[A-Za-z0-9_]{20,}|AKIA[0-9A-Z]{16}|BEGIN (RSA |OPENSSH |EC )?PRIVATE KEY|password[[:space:]]*=[[:space:]]*[^[:space:]]{6,}|authorization:[[:space:]]*(basic|bearer)[[:space:]]+[A-Za-z0-9])" . --glob "!.git/**" --glob "!tests/validate_repository.sh" || true)
if [[ -z "${secret_matches}" ]]; then
  record "Credential patterns" "PASS" "No obvious credential or private-key pattern was found" "None" "None"
else
  redacted_files=$(printf '%s\n' "${secret_matches}" | tr '\n' ' ')
  record "Credential patterns" "FAIL" "Potential secret pattern found; values redacted" "None" "${redacted_files}"
fi

raw_files=$(find . -path "./.git" -prune -o -type f \( -iname "*.rds" -o -iname "*.rdata" -o -iname "*.rda" -o -iname "*.fastq" -o -iname "*.fastq.gz" -o -iname "*.fq" -o -iname "*.fq.gz" -o -iname "*.bam" -o -iname "*.h5" -o -iname "*.h5ad" \) -print)
if [[ -z "${raw_files}" ]]; then
  record "Raw-data files" "PASS" "No excluded research-data extension is present" "None" "None"
else
  record "Raw-data files" "FAIL" "Research-data files are present" "None" "${raw_files//$'\n'/ }"
fi

large_files=$(find . -path './.git' -prune -o -type f -size +99M -print)
if [[ -z "${large_files}" ]]; then
  record "GitHub file-size limit" "PASS" "No repository file exceeds 99 MiB" "None" "None"
else
  record "GitHub file-size limit" "FAIL" "Oversized files are present" "None" "${large_files//$'\n'/ }"
fi

internal_paths=$(rg -l -i "(/data/(projects|gpfs)/|punim[0-9]{4}|[A-Za-z0-9._%+-]+@unimelb\.edu\.au)" . --glob "!.git/**" --glob "!tests/validate_repository.sh" || true)
if [[ -z "${internal_paths}" ]]; then
  record "Institutional paths and identities" "PASS" "No blocked internal path, allocation, or email pattern was found" "None" "None"
else
  record "Institutional paths and identities" "FAIL" "Public files contain blocked internal information" "None" "${internal_paths//$'\n'/ }"
fi

non_english=$(rg -l -P "[\p{Han}]" . --glob "*.md" --glob "*.R" --glob "*.sh" --glob "*.yml" --glob "*.tsv" --glob "!.git/**" || true)
if [[ -z "${non_english}" ]]; then
  record "English public text" "PASS" "No Han-script character was found in public text/code" "None" "This is a character-level check, not editorial certification"
else
  record "English public text" "FAIL" "Non-English characters were found" "None" "${non_english//$'\n'/ }"
fi

non_executable=$(find scripts tests -type f \( -name '*.sh' -o -name '*.R' \) ! -perm -u+x -print)
if [[ -z "${non_executable}" ]]; then
  record "Executable permissions" "PASS" "All Shell and R scripts are user-executable" "None" "None"
else
  record "Executable permissions" "FAIL" "Executable permission is missing" "None" "${non_executable//$'\n'/ }"
fi

gitignore_missing=()
for pattern in '*.rds' '*.RData' '*.fastq.gz' '*.bam' 'outputs/' 'work/' '*.log' '.Rhistory' '*.pem'; do
  grep -Fqx "${pattern}" .gitignore || gitignore_missing+=("${pattern}")
done
if (( ${#gitignore_missing[@]} == 0 )); then
  record ".gitignore coverage" "PASS" "Required data, output, log, temporary, R-state, and credential patterns are excluded" "None" "None"
else
  record ".gitignore coverage" "FAIL" "Required ignore patterns are absent" "None" "${gitignore_missing[*]}"
fi

manifest_missing=()
while IFS=$'\t' read -r figure_file public_path _rest; do
  [[ "${figure_file}" == "figure_file" ]] && continue
  [[ -z "${public_path}" || "${public_path}" == "NA" ]] && continue
  [[ -e "${public_path}" ]] || manifest_missing+=("${public_path}")
done < docs/FIGURE_MANIFEST.tsv
if (( ${#manifest_missing[@]} == 0 )); then
  record "Figure manifest paths" "PASS" "All 42 archived figure paths exist" "None" "None"
else
  record "Figure manifest paths" "FAIL" "Manifest paths are missing" "None" "${manifest_missing[*]}"
fi

if sha256sum -c figures/SHA256SUMS >/dev/null 2>&1; then
  record "Figure checksums" "PASS" "All 42 archived figure checksums match figures/SHA256SUMS" "None" "None"
else
  record "Figure checksums" "FAIL" "At least one archived figure checksum does not match" "None" "Restore the verified figure file or checksum"
fi

if command -v python3 >/dev/null 2>&1 &&
  python3 -c 'import yaml' >/dev/null 2>&1; then
  if python3 -c 'import yaml; yaml.safe_load(open("environment/environment.yml")); yaml.safe_load(open("CITATION.cff"))'; then
    record "YAML syntax" "PASS" "environment.yml and CITATION.cff parse with PyYAML" "None" "None"
  else
    record "YAML syntax" "FAIL" "A YAML file failed to parse" "None" "Correct YAML syntax"
  fi
else
  record "YAML syntax" "NOT TESTED" "Python with PyYAML is unavailable" "None" "Validate with a YAML parser"
fi

if (( failures > 0 )); then
  record "Overall validation" "FAIL" "${failures} failed check(s), ${warnings} warning check(s)" "None" "Resolve all FAIL rows before publication"
  exit 1
fi

record "Overall validation" "PASS" "No failed checks; ${warnings} warning check(s)" "None" "Review NOT TESTED and WARNING rows"
