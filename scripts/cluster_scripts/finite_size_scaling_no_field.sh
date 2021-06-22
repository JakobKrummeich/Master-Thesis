#!/bin/bash 
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/simulation/semi_gcmc_no_field

srun --ntasks=1 --error=error_stream_output/N=16000_%J.err ./semi_gcmc_16000 0.6 1.0 0.02 0.5 300000 0.1 20 20 &

wait
