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
number_of_equilibration_sweeps=$3
number_of_data_taking_sweeps=$4
totalNumberOfParticles=$5

datapath="data/MDSimulations/equilibrate_SGC_MD/N=${totalNumberOfParticles}/T=${temperature}" # relative to root directory of repository

new_directory="../../../${datapath}"
make_directory_if_necessary new_directory
new_directory="${new_directory}/singleRunData"
make_directory_if_necessary new_directory

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=normal,oip
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/MDSimulation/equilibrate_SGC_MD"


for (( runNumber=0; runNumber<${number_of_runs}; runNumber++ )); do

	result_directory="${new_directory}/${runNumber}/"
	make_directory_if_necessary result_directory

	if [ -d "${result_directory}" ]; then
		if [ ! "$(ls -A ${result_directory})" ]; then

			submit_filename="equilibrate_SGC_MD_N=${totalNumberOfParticles}_T=${temperature}_${runNumber}.sh"
			echo "$settings" > ${submit_filename}
			echo -e >> ${submit_filename}

			srun_command="srun --ntasks=1 ./equilibrate_SGC_MD_N=${totalNumberOfParticles} ${temperature} ${result_directory} ${number_of_equilibration_sweeps} ${number_of_data_taking_sweeps} &

		wait"
			echo "$srun_command" >> ${submit_filename}
			sbatch ${submit_filename}
			rm ${submit_filename}

		fi
	fi

done

