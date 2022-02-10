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

for temperature_dir in ${sourcepath}/T=*; do

	for runDir in ${temperature_dir}/*; do

		velocityFiles+=(${runDir}/avgVelocities*.dat)
	
	done

	[[ ${temperature_dir} =~ /T=([[:digit:]]*.[[:digit:]]*) ]] && temperature=${BASH_REMATCH[1]}

	echo 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2
	echo "T=${temperature} analysis running." 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2

	filename=${velocityFiles[0]}

	[[ ${filename} =~ N=([[:digit:]]*) ]] && number_of_particles=${BASH_REMATCH[1]}

	results=$(echo ${velocityFiles[@]} | ./combine_velocities ${target_filepath})

done
