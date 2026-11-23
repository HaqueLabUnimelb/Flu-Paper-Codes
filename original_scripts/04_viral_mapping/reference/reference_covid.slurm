#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Wed Oct 23 2024 09:20:22 GMT+1100 (Australian Eastern Daylight Time)

#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --job-name="spleen_stomics_SAW"
#SBATCH --account="punim2394"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --mail-user=<email>
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --time=31-8:00:00


cd /data/projects/punim2394/Virus_ref_genome/generate_seeker_reference

module load Anaconda3/2023.07-2
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate star_2.6.1d

echo "[DEBUG] CONDA_PREFIX=$CONDA_PREFIX"
hash -r || true
echo "[DEBUG] PATH=$PATH"
echo "[DEBUG] command -v STAR=$(command -v STAR)"
type -a STAR || true
STAR --version || true


bash generate_seeker_reference.sh \
  /data/projects/punim2394/Virus_ref_genome/generate_seeker_reference/combined_GRCm38_SARS2/GRCm38_plus_SARSCoV2.fa \
  /data/projects/punim2394/Virus_ref_genome/generate_seeker_reference/combined_GRCm38_SARS2/GRCm38_plus_SARSCoV2.gtf \
  Mt \
  GRCm38_covid19
