#!/bin/bash

sourcepath=$1

[[ ${sourcepath} =~ (N=[[:digit:]]*) ]] && number_of_particles_string=${BASH_REMATCH[1]}
phase_diagram_filename="${sourcepath}/analyzed_data/phase_diagram_${number_of_particles_string}.dat"

printf "T\tfirst_moment_xA\n" > ${phase_diagram_filename}

for temperature_dir in ${sourcepath}/T=*; do

	NA_files=(${temperature_dir}/NA*.dat)
	
	[[ ${temperature_dir} =~ T=([[:digit:]]*.[[:digit:]]*) ]] && temperature=${BASH_REMATCH[1]}

	filename=${NA_files[0]}

	[[ ${filename} =~ NA_Series_(.*_epsAB=[[:digit:]]*.[[:digit:]]*)_ ]] && distribution_filename=${BASH_REMATCH[1]}
	[[ ${filename} =~ N=([[:digit:]]*)_ ]] && number_of_particles=${BASH_REMATCH[1]}


	for NA_file in ${NA_files[@]}; do
		file_pairs+=($NA_file)
		[[ ${NA_file} =~ NA_Series_(.*_epsAB.*) ]] && file_info=${BASH_REMATCH[1]}
		file_pairs+=($sourcepath/PotEnergySeries_${file_info})
	done

	first_moment=$(echo ${file_pairs[@]} | ./compute_distribution_and_moments $distribution_filename $number_of_particles "${sourcepath}/analyzed_data/")

	printf "${temperature}\t${first_moment}\n" >> ${phase_diagram_filename}

done
