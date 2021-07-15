#!/bin/bash 
#SBATCH --ntasks=1
#SBATCH --partition=oip
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/data_analysis/

srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=16000_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/Master-Thesis/data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh=0.75_long_runs/N=16000 1000000 &

wait
