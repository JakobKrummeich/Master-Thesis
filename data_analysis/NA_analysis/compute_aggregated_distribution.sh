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

[[ ${sourcepath} =~ (N=[[:digit:]]*) ]] && number_of_particles_string=${BASH_REMATCH[1]}
target_filepath="${sourcepath}/aggregated_data/"
make_directory_if_necessary target_filepath
phase_diagram_filename="${target_filepath}first_moment_binder_cumulant_${number_of_particles_string}.dat"

printf "first_moment_xA\tBinder_cumulant\n" > ${phase_diagram_filename}

for runDir in ${sourcepath}/singleRunData/*; do

	NA_files+=(${runDir}/histogram_*.dat)
	
done

filename=${NA_files[0]}

[[ ${filename} =~ N=([[:digit:]]*) ]] && number_of_particles=${BASH_REMATCH[1]}
distribution_filename="combined_N=${number_of_particles}"

results=$(echo ${NA_files[@]} | ./combine_histograms $distribution_filename $number_of_particles ${target_filepath})

printf "${results}\n" >> ${phase_diagram_filename}

