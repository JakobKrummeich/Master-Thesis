#!/bin/bash

sourcepath=$1

NA_files=($sourcepath/NA*.dat)

filename=${NA_files[0]}

[[ ${filename} =~ NA_Series_(.*_epsAB=[[:digit:]]*.[[:digit:]]*)_ ]] && extracted_file_info=${BASH_REMATCH[1]}
[[ ${filename} =~ N=([[:digit:]]*)_ ]] && number_of_particles=${BASH_REMATCH[1]}


for NA_file in ${NA_files[@]}
do
	file_pairs+=($NA_file)
	[[ ${NA_file} =~ NA_Series_(.*_epsAB.*) ]] && file_info=${BASH_REMATCH[1]}
	file_pairs+=($sourcepath/PotEnergySeries_${file_info})
done

echo ${file_pairs[@]} | ./analyze_series $extracted_file_info $number_of_particles

