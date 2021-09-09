#!/bin/bash

sourcepath=$1
min_number_of_equilibration=$2
restrict_to_half_dist=$3

function make_directory_if_necessary() {
	echo -n "${!1} "
	if [ -d "${!1}" ]; then
		echo "exists."
	else
		echo "does not exist. Creating it."
		mkdir ${!1}
	fi
}

[[ ${sourcepath} =~ (N=[[:digit:]]*) ]] && number_of_particles_string=${BASH_REMATCH[1]}
target_filepath="${sourcepath}/analyzed_data/"
make_directory_if_necessary target_filepath
phase_diagram_filename="${target_filepath}structure_factor_from_fluctuations_${number_of_particles_string}.dat"

printf "T\tScc(0)\n" > ${phase_diagram_filename}

for temperature_dir in ${sourcepath}/T=*; do

	NA_files=(${temperature_dir}/NA*.dat)

	[[ ${temperature_dir} =~ /T=([[:digit:]]*.[[:digit:]]*) ]] && temperature=${BASH_REMATCH[1]}

	echo 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2
	echo "T=${temperature} analysis running." 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2

	filename=${NA_files[0]}

	[[ ${filename} =~ NA_Series_(.*_epsAB=[[:digit:]]*.[[:digit:]]*)_ ]] && distribution_filename=${BASH_REMATCH[1]}
	[[ ${filename} =~ N=([[:digit:]]*)_ ]] && number_of_particles=${BASH_REMATCH[1]}

	results=$(echo ${NA_files[@]} | ./compute_structure_factor_from_fluctuations $distribution_filename $number_of_particles ${target_filepath} ${min_number_of_equilibration} ${restrict_to_half_dist})

	printf "${temperature}\t${results}\n" >> ${phase_diagram_filename}

done
