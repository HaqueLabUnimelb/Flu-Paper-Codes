#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Mon Mar 25 2024 15:41:42 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=cascade,sapphire

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="lung_curioseeker_qc"

# The project ID which this job should run under:
#SBATCH --account="punim1236"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The amount of memory in megabytes per node:
#SBATCH --mem=8192

# Use this email address:
#SBATCH --mail-user=<email>

# Send yourself an email when the job:
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=14-0:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from this directory:
source ~/.bashrc
cd /data/scratch/projects/punim1236/curioseeker-2.0.0

# The job command(s):
conda_load
conda activate singularity

## run nf
nextflow run main.nf \
        --input /data/scratch/projects/punim1236/lung_curio/samplesheet_igenome_lung.csv \
        --outdir /data/scratch/projects/punim1236/lung_curio/results_opt_mouse_pr8 \
        -work-dir /data/scratch/projects/punim1236/lung_curio/temp \
        --igenomes_base /data/scratch/projects/punim1236/curioseeker-2.0.0/references/ \
	-resume \
       -profile singularity \
        -config /data/scratch/projects/punim1236/lung_curio/curioseeker_slurm_lung.config

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s
