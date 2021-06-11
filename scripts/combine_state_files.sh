#!/bin/bash

sourcepath=$1

for temperature_dir in ${sourcepath}/T=*; do

	echo -n "Temperature directory: "
	echo ${temperature_dir}

	files=(${temperature_dir}/State*.dat)

	filename=${files[0]}

	[[ ${filename} =~ State_(.*_epsAB=[[:digit:]]*.[[:digit:]]*)_ ]] && condensed_filename="States_${BASH_REMATCH[1]}.dat"

	for file in ${files[@]}; do
		cat $file >> ${temperature_dir}/${condensed_filename}
		rm $file
	done

done
