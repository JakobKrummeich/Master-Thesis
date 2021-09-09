settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=phi
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/data_analysis/"

for N in {16000,8000,4000,2000,1000,500}; do

	submit_filename="compute_phase_diagram_N=${N}.sh"
	echo "$settings" > ${submit_filename}

	srun_command="srun --ntasks=1 --error=error_stream_output/phase_diagram_computation_N=${N}_%J.err ./compute_distributions_and_phase_diagram.sh /home1/krummeich/Master-Thesis/data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh=0.75_long_runs_2/N=${N} 0 &

wait"

	echo "${srun_command}" >> ${submit_filename}
	#sbatch ${submit_filename}
done
