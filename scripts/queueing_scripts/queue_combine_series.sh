settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/scripts/helper_scripts"


for N in {16000,8000,4000,2000}; do
	submit_filename="combine_series_N=${N}.sh"
	echo "$settings" > ${submit_filename}
	echo  -e  >> ${submit_filename}
	srun_command="srun --ntasks=1 --error=error_stream_output/N=${N}_%J.err ./combine_series_files.sh ../../data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh\=0.75_long_runs_3/N\=${N}/ ../../data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh\=0.75_long_runs_2/N\=${N}/  &

wait"
	echo "$srun_command" >> ${submit_filename}
	sbatch ${submit_filename}
done
