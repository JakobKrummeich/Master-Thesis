sourcepath=$1

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/data_analysis/structure_factor_analysis/"


for N in {16000,8000,4000,2000,1000,500}; do

	submit_filename="structure_factors_from_fluctuations_N=${N}.sh"
	echo "$settings" > ${submit_filename}

	srun_command="srun --ntasks=1 --error=error_stream_output/structure_factors_from_fluctuations_%J.err ./compute_structure_factors_from_fluctuations.sh ${sourcepath}/N=${N} 0 0.7469 &

	wait"

	echo "${srun_command}" >> ${submit_filename}
	#sbatch ${submit_filename}

done
