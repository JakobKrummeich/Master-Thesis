#!/bin/bash

sourcepath=$1

for number_of_particles_dir in ${sourcepath}/N=*; do

	for temperature_dir in ${number_of_particles_dir}/T=*; do

		echo -n "Temperature directory: "
		echo ${temperature_dir}

		files=(${temperature_dir}/States_*.dat)

		filepath=${files[0]}

		[[ ${filepath} =~ (States_.*_epsAB=[[:digit:]]*.[[:digit:]]*).dat ]] && filename="${BASH_REMATCH[1]}"

		new_filepath="${temperature_dir}/${filename}_cut.dat"

		echo "${filepath}"
		echo "${new_filepath}"

		cut -f2- "${filepath}" > "${new_filepath}"

		rm "${filepath}"

		mv "${new_filepath}" "${filepath}"

	done

done
