sourcepath=$1

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/data_analysis/"


for temperature_dir in ${sourcepath}/T=*; do
	[[ ${temperature_dir} =~ /T=([[:digit:]]*.[[:digit:]]*) ]] && temperature=${BASH_REMATCH[1]}
	filename="StructureFactor_T=${temperature}.sh"
	echo "$settings" > $filename
	echo  -e  >> $filename

	state_filename=($temperature_dir/States*.dat)
	[[ ${state_filename} =~ States_(.*_epsAB=[[:digit:]]*.[[:digit:]]*).dat ]] && extracted_file_info=${BASH_REMATCH[1]}

	srun_command="srun --ntasks=1 --error=error_stream_output/StructureFactor_computation_T=${temperature}_%J.err ./structure_factor ${state_filename} ${extracted_file_info} 1.0 25.0 16000 200 &

wait"

	echo "$srun_command" >> $filename
	sbatch $filename
done

