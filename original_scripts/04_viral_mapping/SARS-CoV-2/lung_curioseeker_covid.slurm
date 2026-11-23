#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Wed Oct 23 2024 09:20:22 GMT+1100 (Australian Eastern Daylight Time)

#SBATCH -p long
#SBATCH --nodes=1
#SBATCH --job-name="curioseeker"
#SBATCH --account="punim2394"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --mail-user=<email>
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --time=31-00:00:00


# Run the job from this directory:
source ~/.bashrc
cd /data/projects/punim2394/curioseeker_v3.1.0/seeker-v3.1.0

module load GCCcore/13.3.0
module load Nextflow/23.04.2
module load Apptainer/1.4.4

## run nf
nextflow run /data/projects/punim2394/curioseeker_v3.1.0/seeker-v3.1.0/main.nf \
        --input   /data/projects/punim2394/curioseeker_v3.1.0/files_HKX31_covid.csv \
        --outdir /data/projects/punim2394/curioseeker_v3.1.0/result_day4_10_covid\
        -work-dir /data/projects/punim2394/curioseeker_v3.1.0/temp_3 \
        --igenomes_base /data/projects/punim2394/Virus_ref_genome/generate_seeker_reference \
        -profile singularity \
        -config /data/projects/punim2394/curioseeker_v3.1.0/Bioinfo_Seeker_example_slurm.config

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s