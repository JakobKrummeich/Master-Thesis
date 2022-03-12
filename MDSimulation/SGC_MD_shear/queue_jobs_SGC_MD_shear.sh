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
numberOfParticles=$7


targetpath="data/MDSimulations/SGC_MD_shear/N=${numberOfParticles}/T=${temperature}" # relative to root directory of repository

newDirectory="../../../${targetpath}"
make_directory_if_necessary newDirectory
newDirectory+="/shearRate=${shearRate}"
make_directory_if_necessary newDirectory
newDirectory+="/singleRunData"
make_directory_if_necessary newDirectory


for dir in ${initStateDirectory}*/ ; do

	[[ ${dir} =~ T=[[:digit:]]*.[[:digit:]]*/singleRunData/([[:digit:]]*)/ ]] && runNumber="${BASH_REMATCH[1]}"

	result_directory="${newDirectory}/${runNumber}/"
	make_directory_if_necessary result_directory

done

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=normal,oip
#SBATCH --nodes=1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/MDSimulation/SGC_MD_shear"


for dir in ${initStateDirectory}*/ ; do

	[[ ${dir} =~ T=[[:digit:]]*.[[:digit:]]*/singleRunData/([[:digit:]]*)/ ]] && runNumber="${BASH_REMATCH[1]}"

	result_directory="${newDirectory}/${runNumber}/"

	initialStateFile="${dir}N=${numberOfParticles}_final_state.dat"

	srun_command="srun --ntasks=1 ./SGC_MD_shear_N=${numberOfParticles} ${temperature} ${shearRate} ${initialStateFile} ${result_directory} ${numberOfEquilibrationSweeps} ${numberOfDataTakingSweeps} ${stepsize} &

wait"

	submit_filename="SGC_MD_shear_N=${numberOfParticles}_T=${temperature}_shearRate=${shearRate}_${runNumber}.sh"
	echo "$settings" > ${submit_filename}
	echo $'\n' >> ${submit_filename}
	echo "$srun_command" >> ${submit_filename}

	sbatch ${submit_filename}
	rm ${submit_filename}

	sleep 0.02

done


