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

printf "T\tfirst_moment_xA\tBinder_cumulant\n" > ${phase_diagram_filename}

for temperature_dir in ${sourcepath}/T=*; do

	for runDir in ${temperature_dir}/*; do

		NA_files+=(${runDir}/histogram_*.dat)
	
	done

	[[ ${temperature_dir} =~ /T=([[:digit:]]*.[[:digit:]]*) ]] && temperature=${BASH_REMATCH[1]}

	echo 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2
	echo "T=${temperature} analysis running." 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2

	filename=${NA_files[0]}

	[[ ${filename} =~ N=([[:digit:]]*) ]] && number_of_particles=${BASH_REMATCH[1]}
	distribution_filename="combined_N=${number_of_particles}_T=${temperature}"

	results=$(echo ${NA_files[@]} | ./combine_histograms $distribution_filename $number_of_particles ${target_filepath})

	printf "${temperature}\t${results}\n" >> ${phase_diagram_filename}

done
