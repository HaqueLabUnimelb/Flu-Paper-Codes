# Flu Paper Codes

**AI-assisted reproduction:** [Open the Codex reproduction prompt](PROMPTS/CODEX_REPRODUCTION_PROMPT.md)

## Overview

This repository archives analysis scripts associated with a flu manuscript
revision. It provides code for Curio Seeker preprocessing, spatial feature
plots, RCTD cell-type mapping, RCTD visualization, and viral-reference mapping.
The repository does not contain raw data, processed RDS objects, reference
genomes, Curio Seeker itself, or generated manuscript figures.

Full execution requires the input datasets and reference resources described in
[Input data requirements](docs/INPUT_DATA_REQUIREMENTS.md). The scripts were
syntax-checked but were not rerun in full because the workflows require
project data, external references, specialized software, and HPC resources.
This repository should therefore be described as a code archive with
reproduction guidance, not as a fully self-contained reproducibility package.

## Workflow

1. Run Curio Seeker preprocessing for the spatial-transcriptomics samples.
2. Generate spatial quality-control and gene-expression feature plots.
3. Run day-4 and day-10 RCTD as independent large-memory jobs.
4. Inspect an RCTD result with the visualization script.
5. Run HKx31 or SARS-CoV-2 reference mapping when the corresponding
   samplesheets and references are available.

See [Analysis workflow](docs/ANALYSIS_WORKFLOW.md) for details.

## Repository structure

| Path | Contents |
|---|---|
| `scripts/01_preprocessing/` | Curio Seeker preprocessing SLURM script |
| `scripts/02_feature_plots/` | Spatial feature-plot R script |
| `scripts/03_rctd/` | Day-4/day-10 RCTD scripts and visualization |
| `scripts/04_viral_mapping/` | HKx31 and SARS-CoV-2 mapping scripts |
| `config/` | Example path and sample configuration |
| `environment/` | Confirmed and inferred dependency information |
| `docs/` | Workflow, provenance, changes, and troubleshooting |
| `data/` | Input-placement guidance; no data are committed |
| `tests/` | Lightweight repository validation |
| `PROMPTS/` | Standalone Codex reproduction prompt |

## Code-to-analysis and output map

| Analysis | Main script | Expected inputs | Expected outputs | Figure mapping |
|---|---|---|---|---|
| Curio Seeker preprocessing | `scripts/01_preprocessing/lung_curioseeker_qc.sh` | Curio Seeker workflow, samplesheet, reference configuration | Curio Seeker output directory | No direct figure output |
| Spatial feature plots | `scripts/02_feature_plots/FeaturePlots_Flu.R` | Four QCed Seurat RDS objects | JPEG QC and expression plots | Manuscript figure numbers unresolved |
| Day-4 RCTD **(300 GB HPC job)** | `scripts/03_rctd/RCTD_flu_d4_multi_spacexr2.2.R` via `mouseDay4RCTD_full_v2.sh` | PR8 reference and two edited day-4 Seurat RDS objects | `Day4mouseRCTDoutmin30_multi.rds` | No direct figure output |
| Day-10 RCTD **(300 GB HPC job)** | `scripts/03_rctd/RCTD_flu_d10_multi_spacexr2.2.R` via `mouseDay10RCTD_full_v3.sh` | PR8 reference and two edited day-10 Seurat RDS objects | `Day10mouseRCTDoutmin30_multi.rds` | No direct figure output |
| RCTD visualization | `scripts/03_rctd/RCTD_flu_vis.R` | Multi-replicate RCTD RDS result | Interactive plot; optional file when configured | Manuscript figure number unresolved |
| HKx31 mapping | `scripts/04_viral_mapping/lung_curioseeker_HKx31genomemapping.sh` | Curio Seeker workflow, HKx31 samplesheet and references | Configured Nextflow output directory | No direct figure output |
| SARS-CoV-2 mapping | `scripts/04_viral_mapping/lung_curioseeker_SARSCOV2genomemapping.sh` | Curio Seeker workflow, SARS-CoV-2 samplesheet and references | Configured Nextflow output directory | No direct figure output |

The static code-to-figure record is in
[CODE_TO_FIGURE_MAP.md](docs/CODE_TO_FIGURE_MAP.md). No manuscript figures are
committed because their manuscript mapping was not established from the code
alone.

## Requirements

- Linux and Bash
- SLURM for the supplied batch scripts
- Nextflow plus a Singularity-compatible runtime for Curio Seeker workflows
- Curio Seeker workflow files and the relevant reference resources
- R 4.3.x, Seurat 5.1.0, and spacexr 2.2.0 for RCTD
- The additional R packages loaded by the plotting scripts

Confirmed, inferred, and unavailable versions are distinguished in
[environment/README.md](environment/README.md).

## Installation

```bash
git clone https://github.com/HaqueLabUnimelb/Flu-Paper-Codes.git
cd Flu-Paper-Codes
conda env create --prefix environment/rctd_spacexr_2.2 --file environment/environment.yml
```

The environment file prepares the confirmed RCTD dependencies but does not
install spacexr itself. Install spacexr 2.2.0 from the verified source described
in [environment/README.md](environment/README.md). Curio Seeker and its
container images are external requirements.

## Configuration

```bash
cp config/paths.example.env config/paths.env
# Edit config/paths.env with local paths; this file is ignored by Git.
set -a
source config/paths.env
set +a
```

Do not edit scripts to insert institutional paths. Required filenames and
sample labels are summarized in `config/samples.example.tsv`.

## Example execution

Create the log directory before submitting SLURM scripts:

```bash
mkdir -p logs
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/01_preprocessing/lung_curioseeker_qc.sh

sbatch --account=YOUR_ACCOUNT --partition=long \
  scripts/03_rctd/mouseDay4RCTD_full_v2.sh
sbatch --account=YOUR_ACCOUNT --partition=long \
  scripts/03_rctd/mouseDay10RCTD_full_v3.sh

Rscript --vanilla scripts/02_feature_plots/FeaturePlots_Flu.R
Rscript --vanilla scripts/03_rctd/RCTD_flu_vis.R
```

Use site-appropriate accounts and partitions. Never copy an account value from
another institution. The RCTD scripts request 300 GB RAM and should not be run
on a login node. More detail is available in
[HPC execution](docs/HPC_EXECUTION.md).

## Expected outputs

- Curio Seeker and viral mapping: pipeline-specific files under the configured
  output directories.
- Feature plotting: JPEG files under `FLU_FEATURE_PLOT_OUTPUT_DIR`.
- RCTD: one multi-replicate RDS file per time point under
  `FLU_RCTD_OUTPUT_DIR`.
- RCTD visualization: an on-screen plot by default, or a file when
  `FLU_RCTD_VIS_OUTPUT` is set.

## Validation and limitations

Run:

```bash
FLU_RSCRIPT=/path/to/Rscript tests/validate_repository.sh
```

Validation covers repository structure, links, syntax, obvious secret patterns,
internal paths, raw-data extensions, permissions, and GitHub file-size limits.
It does not validate biological interpretation or reproduce analytical results.
See [Reproducibility](docs/REPRODUCIBILITY.md) and the latest
[validation summary](logs/validation_summary.txt).

## Data and code availability

No raw or processed research data are distributed in this repository. Data
availability and access conditions must be supplied by the manuscript authors;
see [data/README.md](data/README.md). The public code is available from this
GitHub repository under the MIT License.

## Citation

Use [CITATION.cff](CITATION.cff) to cite this software archive. The manuscript
citation was not available during repository preparation and should be added by
the authors when confirmed.

## License and contact

Code is released under the [MIT License](LICENSE). For questions or reproducible
bug reports, open an issue in this GitHub repository. Requests involving data
access should follow the manuscript's data-availability statement.
