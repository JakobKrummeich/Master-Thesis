function make_directory_if_necessary() {
	echo -n "${!1} "
	if [ -d "${!1}" ]; then
		echo "exists."
	else 
		echo "does not exist. Creating it."
		mkdir ${!1}
	fi
}

temperature=$1
initStateDirectory=$2

targetpath="data/MDSimulations/SGC_MD_shear/N=1000/T=${temperature}" # relative to root directory of repository

new_directory="../../../${targetpath}"

make_directory_if_necessary new_directory

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/MDSimulation/SGC_MD_shear"


for dir in ${initStateDirectory}*/ ; do

	[[ ${dir} =~ T=[[:digit:]]*.[[:digit:]]*/([[:digit:]]*)/ ]] && runNumber="${BASH_REMATCH[1]}"

	submit_filename="SGC_MD_shear_N=1000_T=${temperature}_${runNumber}.sh"
	echo "$settings" > ${submit_filename}
	echo -e >> ${submit_filename}

	result_directory="${new_directory}/${runNumber}/"
	make_directory_if_necessary result_directory

	initialStateFile="${dir}final_state.dat"

	srun_command="srun --ntasks=1 --error=${result_directory}error_stream.err ./SGC_MD_shear ${temperature} ${initialStateFile} ${result_directory} &

wait"
	echo "$srun_command" >> ${submit_filename}
	sbatch ${submit_filename}
	rm ${submit_filename}

done

