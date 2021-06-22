#!/bin/bash 
#SBATCH --ntasks=6
#SBATCH --partition=oip
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/data_analysis/

srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=500_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/simulation/semi_gcmc_no_field/sorted_data/finite_size_scaling_no_field/epsAB=0.100000/Roh=0.600000/N=500 10000 &

srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=1000_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/simulation/semi_gcmc_no_field/sorted_data/finite_size_scaling_no_field/epsAB=0.100000/Roh=0.600000/N=1000 10000 &

srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=2000_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/simulation/semi_gcmc_no_field/sorted_data/finite_size_scaling_no_field/epsAB=0.100000/Roh=0.600000/N=2000 10000 &

srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=4000_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/simulation/semi_gcmc_no_field/sorted_data/finite_size_scaling_no_field/epsAB=0.100000/Roh=0.600000/N=4000 10000 &

srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=8000_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/simulation/semi_gcmc_no_field/sorted_data/finite_size_scaling_no_field/epsAB=0.100000/Roh=0.600000/N=8000 10000 &

srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=16000_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/simulation/semi_gcmc_no_field/sorted_data/finite_size_scaling_no_field/epsAB=0.100000/Roh=0.600000/N=16000 10000 &

wait
