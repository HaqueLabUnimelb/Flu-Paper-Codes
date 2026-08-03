# Troubleshooting

## A required variable is not set

Copy `config/paths.example.env` to `config/paths.env`, replace all
placeholders, and submit from the repository root. Set `FLU_CONFIG` if the
configuration file is elsewhere.

## A required RDS file is missing

Check the filenames in `config/samples.example.tsv`. Use the explicit RCTD
variables when edited inputs have different names. Do not rename or reassign
samples without author confirmation.

## R or package version checks fail

The RCTD wrappers require R 4.3.x, Seurat 5.1.0, and spacexr 2.2.0. Point
`FLU_R_ENV_PREFIX` to the correct environment. Do not bypass these checks
without documenting a version change.

## tidyverse or mltools is missing

The inspected RCTD environment did not contain these plotting dependencies.
Install them in a separate authorized environment or extend a copy of the
environment, recording the installed versions.

## Nextflow or a container runtime is unavailable

Load site-provided modules or activate an authorized environment before
submitting the Curio Seeker scripts. Confirm that the configured Curio Seeker
release supports the supplied samplesheet and `--igenomes_base` value.

## SLURM cannot open the log file

Run `mkdir -p logs` in the repository root before `sbatch`. SLURM resolves
the log path before the job body starts.

## RCTD requires too much memory

Use a large-memory compute partition. Do not test the full workflow on a login
node. A lightweight input compatibility path is available by setting
`RCTD_PREFLIGHT_ONLY=1`, but it still reads the RDS inputs and may require a
compute node.
