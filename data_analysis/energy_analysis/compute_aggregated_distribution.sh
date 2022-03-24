#!/bin/bash

function make_directory_if_necessary() {
	echo -n "${!1} "
	if [ -d "${!1}" ]; then
		echo "exists."
	else
		echo "does not exist. Creating it."
		mkdir ${!1}
	fi
}

sourcepath=$1
numberOfBins=$2
numberOfEquilibrationLines=$3

[[ ${sourcepath} =~ (N=[[:digit:]]*) ]] && number_of_particles_string=${BASH_REMATCH[1]}
[[ ${number_of_particles_string} =~ N=([[:digit:]]*) ]] && number_of_particles=${BASH_REMATCH[1]}

target_filepath="${sourcepath}/aggregated_data/"
make_directory_if_necessary target_filepath

momentResultFile="${target_filepath}momentsOfAggregatedEnergyDistributions_${number_of_particles_string}.dat"

printf "energyName\tmean\tstandardDeviation\n" > ${momentResultFile}

for runDir in ${sourcepath}/singleRunData/*; do

	energyFiles+=(${runDir}/energySeries.dat)
	
done

results=$(echo ${energyFiles[@]} | ./compute_aggregated_energy_prob_dist ${numberOfBins} ${target_filepath} ${numberOfEquilibrationLines})

printf "${results}" >> ${momentResultFile}

