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
targetpath=$8
maxNumber=$9

targetpath="${targetpath}/N=${numberOfParticles}/T=${temperature}" # relative to root directory of repository

newDirectory="${targetpath}"
make_directory_if_necessary newDirectory
newDirectory+="/shearRate=${shearRate}"
make_directory_if_necessary newDirectory
newDirectory+="/singleRunData"
make_directory_if_necessary newDirectory

count=0
for dir in ${initStateDirectory}*/ ; do

	[[ ${dir} =~ /singleRunData/([[:digit:]]*)/ ]] && runNumber="${BASH_REMATCH[1]}"

	result_directory="${newDirectory}/${runNumber}/"
	make_directory_if_necessary result_directory

	count+=1
	if [ "$count" -ge "$maxNumber" ]; then
		break
	fi

done

count=0
for dir in ${initStateDirectory}*/ ; do

	[[ ${dir} =~ /singleRunData/([[:digit:]]*)/ ]] && runNumber="${BASH_REMATCH[1]}"

	result_directory="${newDirectory}/${runNumber}/"

	if [ -d "${result_directory}" ]; then
		if [ ! "$(ls -A ${result_directory})" ]; then

			initialStateFile="${dir}N=${numberOfParticles}_final_state.dat"

			srun_command="srun --ntasks=1 ./SGC_MD_shear_N=${numberOfParticles} ${temperature} ${shearRate} ${initialStateFile} ${result_directory} ${numberOfEquilibrationSweeps} ${numberOfDataTakingSweeps} ${stepsize} &

wait"

			settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=normal,oip
#SBATCH --nodes=1
#SBATCH --output=${result_directory}/slurm-%j.out
#SBATCH --error=${result_directory}/slurm-%j.err
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --exclude=broadwell09,broadwell47
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/MDSimulation/SGC_MD_shear"

			submit_filename="SGC_MD_shear_N=${numberOfParticles}_T=${temperature}_shearRate=${shearRate}_${runNumber}.sh"
			echo "$settings" > ${submit_filename}
			echo $'\n' >> ${submit_filename}
			echo "$srun_command" >> ${submit_filename}

			sbatch ${submit_filename}
			rm ${submit_filename}

		fi
	fi

	count+=1
	if [ "$count" -ge "$maxNumber" ]; then
		break
	fi

done

