#!/bin/bash

sourcepath=$1
targetpath=$2
included_pattern=$3

function make_directory_if_necessary() {
	echo -n "${!1} "
	if [ -d "${!1}" ]; then
		echo "exists."
	else 
		echo "does not exist. Creating it."
		mkdir ${!1}
	fi
}


for filename in $sourcepath/*${included_pattern}*.dat; do

		[[ ${filename} =~ _epsAB=([[:digit:]]*.[[:digit:]]*). ]] && epsAB=${BASH_REMATCH[1]}
		currentpath=${targetpath}/epsAB=$epsAB
		make_directory_if_necessary currentpath

		[[ ${filename} =~ _AvgDens=([[:digit:]]*.[[:digit:]]*)_ ]] && density=${BASH_REMATCH[1]}
		currentpath=${currentpath}/Roh=$density
		make_directory_if_necessary currentpath

		[[ ${filename} =~ _N=([[:digit:]]*)_ ]] && number_of_particles=${BASH_REMATCH[1]}
		currentpath=${currentpath}/N=$number_of_particles
		make_directory_if_necessary currentpath
		
		[[ ${filename} =~ _T=([[:digit:]]*.[[:digit:]]*)_ ]] && temperature=${BASH_REMATCH[1]}
		currentpath=${currentpath}/T=$temperature
		make_directory_if_necessary currentpath

		mv $filename $currentpath

done
