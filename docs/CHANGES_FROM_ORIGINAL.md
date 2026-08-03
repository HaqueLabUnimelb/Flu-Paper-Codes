# Changes from the supplied scripts

The original files were not modified. All changes below apply only to copied
public versions. No analytical threshold, reference choice, time-point filter,
sample grouping, RCTD mode, or plotting scale was intentionally changed.

| Public file | Original behavior | Public behavior | Reason | Classification |
|---|---|---|---|---|
| `lung_curioseeker_qc.sh` | Used institutional paths, account, partition, email, shell startup file, and local Conda helper; the missing-SLURM branch called `exit1` | Reads paths from a local config, omits identity/allocation directives, checks dependencies, creates outputs, and uses `exit 1` | Safety and portability | Portability/validation |
| HKx31 mapping Shell script | Used institutional paths, account, partition, email, and local Conda helper | Uses configurable paths and external environment preparation; job name identifies HKx31 | Public HPC portability | Portability/cosmetic |
| SARS-CoV-2 mapping Shell script | Used institutional paths, account, partition, email, and local Conda helper | Uses configurable paths and external environment preparation; job name identifies SARS-CoV-2 | Public HPC portability | Portability/cosmetic |
| `FeaturePlots_Flu.R` | Read four absolute RDS paths and wrote to absolute output directories | Uses `FLU_PROCESSED_DATA_DIR` and `FLU_FEATURE_PLOT_OUTPUT_DIR`; checks inputs and creates directories | Portable, informative failure | Portability/validation |
| Day-4 RCTD R script | Read absolute reference/puck paths and wrote to an absolute directory | Uses environment variables with repository-relative defaults; checks inputs/output and coordinate dimensionality | Portable, informative failure | Portability/validation |
| Day-10 RCTD R script | Read absolute reference/puck paths and wrote to an absolute directory | Uses environment variables with repository-relative defaults; checks inputs/output and coordinate dimensionality | Portable, informative failure | Portability/validation |
| Day-4 RCTD Shell wrapper | Used an institutional account, log path, project path, and environment path | Resolves the repository location, reads configuration, uses a relative log path, and allows environment/module overrides | Public HPC portability | Portability |
| Day-10 RCTD Shell wrapper | Used an institutional account, log path, project path, and environment path | Resolves the repository location, reads configuration, uses a relative log path, and allows environment/module overrides | Public HPC portability | Portability |
| `RCTD_flu_vis.R` | Read four results consecutively, so only the final day-10 full-mode assignment was effective; displayed but did not save the plot | Reads one configurable result, defaults to the same effective day-10 output, prints the plot, and optionally saves when explicitly requested | Clarify effective input and enable controlled output | Portability/validation |

Additional repository-wide changes:

- Added English purpose/input/output headers and execution examples.
- Replaced personal and allocation information with user-supplied SLURM options.
- Added quoted variables, `set -euo pipefail`, missing-input checks, and safe
  output-directory creation where applicable.
- Retained the original Nextflow flags, RCTD parameters, feature names, plot
  dimensions, color scales, and explicitly saved output names.
- Archived 42 existing JPEG outputs whose names and locations match explicit
  saves in `FeaturePlots_Flu.R`; the image bytes were not modified.
- No random seed was added because the supplied scripts did not set one and the
  appropriate scientific choice was not established.

The optional RCTD visualization export is operationally new but is disabled by
default. It does not change the plotted data or aesthetics.
