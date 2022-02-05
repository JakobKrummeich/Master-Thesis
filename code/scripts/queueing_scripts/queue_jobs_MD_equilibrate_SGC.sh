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

number_of_runs=$2

datapath="data/MDSimulations/equilibrate_SGC_MD/N=1000/T=${temperature}" # relative to root directory of repository

new_directory="../../../${datapath}"

make_directory_if_necessary new_directory

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/MDSimulation/equilibrate_SGC_MD"


for (( runNumber=0; runNumber<${number_of_runs}; runNumber++ )); do
	submit_filename="equilibrate_SGC_MD_N=1000_T=${temperature}_${runNumber}.sh"
	echo "$settings" > ${submit_filename}
	echo  -e  >> ${submit_filename}

	result_directory="${new_directory}/${runNumber}/"
	make_directory_if_necessary result_directory

	srun_command="srun --ntasks=1 --error=${result_directory}/error_stream.err ./equilibrate_SGC_MD ${temperature} ${result_directory} &

wait"
	echo "$srun_command" >> ${submit_filename}
	sbatch ${submit_filename}
	rm ${submit_filename}

done

