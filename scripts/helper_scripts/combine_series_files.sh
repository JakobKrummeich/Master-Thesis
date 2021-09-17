#!/bin/bash

sourcepath=$1
targetpath=$2

for temperature_dir in ${sourcepath}/T=*; do

	echo -n "Temperature directory: "
	echo ${temperature_dir}

	[[ ${temperature_dir} =~ (T=[[:digit:]]*.[[:digit:]]*) ]] && temperature=${BASH_REMATCH[1]}

	echo ${temperature}

	files=(${temperature_dir}/*Series*.dat)

	for file in ${files[@]}; do
		[[ ${file} =~ T=[[:digit:]]*.[[:digit:]]*/(.*Series_.*.dat) ]] && filename=${BASH_REMATCH[1]}
		echo ${filename}
		echo ${targetpath}/${temperature}/${filename}
		tail -n +2 $file >> ${targetpath}/${temperature}/${filename}
		rm $file
	done

done
