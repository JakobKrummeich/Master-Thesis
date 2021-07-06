#!/bin/bash

sourcepath=$1

for temperature_dir in ${sourcepath}/T=*; do

	echo -n "Temperature directory: "
	echo ${temperature_dir}

	files=(${temperature_dir}/State_*.dat)

	filename=${files[0]}

	[[ ${filename} =~ State_(.*_epsAB=[[:digit:]]*.[[:digit:]]*)_ ]] && condensed_filename="States_${BASH_REMATCH[1]}.dat"

	for file in ${files[@]}; do
		[[ ${file} =~ [[:digit:]]*.[[:digit:]]*_([[:digit:]]_[[:digit:]]).dat ]] && state_label=${BASH_REMATCH[1]}
		echo "State label (first number is RunNumber, second number is number of saved state of that run): ${state_label}"  >> ${temperature_dir}/${condensed_filename}
		cat $file >> ${temperature_dir}/${condensed_filename}
		rm $file
	done

done
