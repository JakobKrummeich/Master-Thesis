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
shearRate=$2
initStateDirectory=$3
numberOfEquilibrationSweeps=$4
numberOfDataTakingSweeps=$5
stepsize=$6


targetpath="data/MDSimulations/SGC_MD_shear/N=1000/T=${temperature}" # relative to root directory of repository

newDirectory="../../../${targetpath}"
make_directory_if_necessary newDirectory
newDirectory+="/shearRate=${shearRate}"
make_directory_if_necessary newDirectory
newDirectory+="/singleRunData"
make_directory_if_necessary newDirectory

numberOfInitialStates=0

for dir in ${initStateDirectory}*/ ; do

	[[ ${dir} =~ T=[[:digit:]]*.[[:digit:]]*/([[:digit:]]*)/ ]] && runNumber="${BASH_REMATCH[1]}"

	result_directory="${newDirectory}/${runNumber}/"
	make_directory_if_necessary result_directory

done

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/MDSimulation/SGC_MD_shear"



for dir in ${initStateDirectory}*/ ; do

	[[ ${dir} =~ T=[[:digit:]]*.[[:digit:]]*/([[:digit:]]*)/ ]] && runNumber="${BASH_REMATCH[1]}"

	result_directory="${newDirectory}/${runNumber}/"

	initialStateFile="${dir}final_state.dat"

	srun_command="srun --ntasks=1 --error=${result_directory}error_stream.err ./SGC_MD_shear ${temperature} ${shearRate} ${initialStateFile} ${result_directory} ${numberOfEquilibrationSweeps} ${numberOfDataTakingSweeps} ${stepsize} &

wait"

	submit_filename="SGC_MD_shear_N=1000_T=${temperature}_shearRate=${shearRate}_${runNumber}.sh"
	echo "$settings" > ${submit_filename}
	echo $'\n' >> ${submit_filename}
	echo "$srun_command" >> ${submit_filename}

	sbatch ${submit_filename}
	rm ${submit_filename}

	sleep 0.05

done


