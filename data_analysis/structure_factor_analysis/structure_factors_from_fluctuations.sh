#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/data_analysis/structure_factor_analysis/
srun --ntasks=1 --error=error_stream_output/structure_factors_from_fluctuations_%J.err ./compute_structure_factors_from_fluctuations.sh ../../data/semi_gcmc_no_field/sorted_data/roh=0.75_missing_temperatures_for_structure_factor_comparison/N=16000/ 0 0 &

wait
