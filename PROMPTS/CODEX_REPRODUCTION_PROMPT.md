# Codex prompt for reproducing the Flu Paper Codes workflow

Use this prompt together with a local clone of the Flu Paper Codes repository.
Fill in the values below before starting. Leave a value blank only when it is
not applicable.

```text
RAW_DATA_DIRECTORY=
OUTPUT_DIRECTORY=
EXECUTION_MODE=SLURM_OR_LOCAL
SLURM_ACCOUNT=
SLURM_PARTITION=
REFERENCE_DIRECTORY=
```

You are responsible for technically reproducing the computational workflow
documented in this repository. Treat the repository documentation and public
scripts as the source of truth. Preserve all documented scientific parameters,
sample assignments, reference choices, thresholds, plotting settings, and RCTD
options. Do not infer or silently change a scientifically consequential value.

## 1. Consolidated technical permission preflight

Before reading data or running commands, ask me once for all technical
permissions required for this session. Combine the request so that the workflow
is not repeatedly interrupted. Request only permissions that are relevant:

- read the repository, raw-data directory, processed-data directory, and
  reference directory;
- identify file formats, filenames, sizes, checksums, and directory structure;
- create the selected output directory and separate work/cache directories;
- run non-destructive Shell and R commands;
- inspect module, Conda, R, package, container, Nextflow, and SLURM information;
- run lightweight syntax, dependency, and input-compatibility checks;
- install missing packages or tools only when necessary and explicitly
  authorized, using an isolated environment;
- read reference genomes and indexes;
- submit and monitor SLURM jobs only when explicitly authorized;
- write command logs, job scripts, checksums, inventories, and validation
  reports under the selected output directory.


If broad but safe approval is available, request it only for the technical
directories and commands needed here. Do not request access to unrelated
projects, home-directory secrets, shell startup files, or raw data outside the
directory I supplied.

## 2. Safety and execution rules

- Never modify, rename, move, or delete raw data.
- Write every generated file to `OUTPUT_DIRECTORY` or an explicitly approved
  work directory.
- Do not upload raw data, processed research objects, credentials, or internal
  configuration to GitHub or another external service.
- Do not run full analyses on an HPC login node.
- Do not submit, cancel, resubmit, or alter SLURM jobs without authorization.
- Avoid high-frequency queue polling.
- Stop safely when a required file is missing.
- Do not guess sample identities, biological groups, time points, references,
  or file substitutions.
- Ask for clarification only when a scientifically consequential decision
  cannot be resolved from repository documentation and supplied metadata.
- Never claim a stage succeeded until its expected files have been verified.

## 3. Inspect the repository and supplied data

Read, in this order:

1. `README.md`
2. `docs/ANALYSIS_WORKFLOW.md`
3. `docs/INPUT_DATA_REQUIREMENTS.md`
4. `docs/HPC_EXECUTION.md`
5. `docs/REPRODUCIBILITY.md`
6. `environment/README.md`
7. `config/paths.example.env`
8. `config/samples.example.tsv`

Inventory the supplied data without changing it. For each relevant file,
record:

- relative path under `RAW_DATA_DIRECTORY`;
- file type and extension;
- size and modification date;
- SHA-256 checksum where practical;
- sample name or identifier explicitly supported by metadata or filename;
- whether it matches a documented required input;
- any ambiguity, duplicate candidate, or missing requirement.

Inspect headers or object metadata only with lightweight, format-appropriate
commands. Do not load every large object at once on a login node. Compare the
inventory with the exact filenames, Seurat metadata fields, SPATIAL reduction,
samplesheets, and references required by the repository. Report mismatches
before running analysis.

## 4. Create local configuration

Copy `config/paths.example.env` to a new configuration file outside the Git
working tree or to the Git-ignored path `config/paths.env`. Populate it with
the supplied directories and verified file assignments. Do not edit public
scripts to insert absolute paths.

Record the final configuration in the reproduction report, but redact private
institutional prefixes if the report will be shared publicly. Never include
tokens, passwords, private keys, or authentication cookies.

## 5. Verify dependencies

Before execution:

- identify the operating system and whether the current host is a login or
  compute node;
- check Bash, SLURM, Nextflow, the Singularity-compatible runtime, Curio Seeker,
  R, and required R packages;
- record exact versions when commands are available;
- compare R, Seurat, and spacexr against the versions enforced by the RCTD
  wrappers;
- verify the configured Curio Seeker workflow and `main.nf`;
- verify samplesheets, Nextflow configuration, references, and indexes;
- run `tests/validate_repository.sh`;
- run `bash -n` for Shell scripts;
- parse R scripts without executing the full workflow.

If a required dependency is missing, stop that stage. Propose an isolated
installation method and wait for authorization before installing. Record both
confirmed and unresolved versions. Do not report an inferred version as
confirmed.

## 6. Choose execution mode

### SLURM

Use the site allocation and partition supplied above, not values from historical
scripts. Create explicit job names and log paths. Follow
`docs/HPC_EXECUTION.md`. The RCTD jobs are documented as one-node,
one-CPU, 300 GB, 31-day large-memory jobs. Create required output and log
directories before submission.

Before each `sbatch):

1. show the exact command and script;
2. verify configured paths;
3. perform lightweight syntax and dependency checks;
4. obtain submission authorization if it was not already granted;
5. record the returned job ID.

Monitor with low-frequency `squeue`, `sacct`, or
`scontrol show job`. Do not treat a submitted or running job as completed.
After completion, inspect the exit code and expected outputs.

### Local Linux

Run only stages that fit available CPU, memory, storage, and runtime limits.
Translate SLURM wrappers into explicit local commands without changing their
analysis arguments. Warn that RCTD may require substantially more memory than a
typical workstation. Do not launch a stage when resources are clearly
insufficient.

## 7. Execute in documented order

Use this order unless a stage is not required or its prerequisite is missing:

1. Curio Seeker preprocessing.
2. Spatial feature plots.
3. Day-4 RCTD.
4. Day-10 RCTD.
5. RCTD visualization.
6. HKx31 and/or SARS-CoV-2 mapping when the corresponding inputs are supplied.

Independent stages may be scheduled separately, but never hide skipped
prerequisites. Log every executed command with timestamp, working directory,
exit status, relevant environment variables, and job ID. Avoid logging secrets.

Do not modify the repository's scientific settings. In particular, retain the
documented time-point filters, sample grouping, UMI minimum, cell-type minimum
instances, doublet mode, feature names, plot dimensions, and color scales.

## 8. Validate after every stage

For each stage:

- confirm the command exit status;
- inspect logs for errors and warnings;
- verify each expected output exists;
- record output size, modification time, and checksum where practical;
- confirm outputs are outside the raw-data directory;
- identify figure files without assigning manuscript panel numbers unless
  supported by repository documentation or author-provided evidence;
- mark the stage `COMPLETED`, `FAILED`, `SKIPPED`, or `NOT TESTED`.

Do not create substitute results for a missing output. Record failed commands
and error messages verbatim where safe, while redacting secrets and private
credentials.

## 9. Final reproduction report

Create a Markdown report under `OUTPUT_DIRECTORY` containing:

1. input inventory;
2. configuration used, with sensitive paths redacted in any public copy;
3. commands executed and exit statuses;
4. software, package, workflow, container, and reference versions;
5. completed stages;
6. failed, skipped, and untested stages;
7. SLURM job IDs and final states, when applicable;
8. generated outputs with sizes and checksums where practical;
9. generated figure files and any evidence-supported mapping;
10. missing inputs, dependencies, or references;
11. deviations from the documented workflow;
12. warnings and unresolved decisions;
13. an explicit statement of whether reproduction was fully verified,
    partially completed, or not completed.

Never claim successful reproduction unless all required stages completed and
their expected outputs were verified. A successful syntax check, job
submission, or file creation is not equivalent to scientific reproduction.
