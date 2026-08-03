# HPC execution

The supplied batch scripts are templates. Review every resource directive and
pass the correct allocation and partition for the target cluster. Never reuse
an account value copied from another institution.

## Preparation

```bash
cp config/paths.example.env config/paths.env
# Edit local paths.
mkdir -p logs
```

Load Nextflow and a Singularity-compatible runtime before submitting Curio
Seeker jobs, or configure site modules in the job environment.

## Example submissions

```bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/01_preprocessing/lung_curioseeker_qc.sh

sbatch --account=YOUR_ACCOUNT --partition=long \
  scripts/03_rctd/mouseDay4RCTD_full_v2.sh

sbatch --account=YOUR_ACCOUNT --partition=long \
  scripts/03_rctd/mouseDay10RCTD_full_v3.sh

sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/04_viral_mapping/lung_curioseeker_HKx31genomemapping.sh

sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/04_viral_mapping/lung_curioseeker_SARSCOV2genomemapping.sh
```

The RCTD wrappers request one node, one CPU, 300 GB RAM, and 31 days in the
`long` partition. They are large-memory jobs. Adjust the partition only when
the target cluster requires it; do not reduce memory without a justified test.

The scripts do not submit, cancel, monitor, or resubmit jobs automatically.
Use low-frequency `squeue`, `sacct`, and `scontrol show job` checks
according to local policy.

## Local Linux

The R scripts can run locally when all inputs, packages, and sufficient memory
are available. The Nextflow launchers and supplied wrappers require SLURM;
translate their commands carefully for a local environment. RCTD may exceed
typical workstation memory.
