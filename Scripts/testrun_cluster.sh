#!/bin/bash
#SBATCH --job-name=semi_gcmc_no_field
#SBATCH --mem=800mb
#SBATCH --partition=test
#SBATCH --nodes=1-1
#SBATCH --ntasks=1
#SBATCH --error=slurm-%A_%a.err

module load intel
./semi_gcmc

