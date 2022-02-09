#!/bin/bash

sourcepath=$1
targetpath=$2

function make_directory_if_necessary() {
	echo -n "${!1} "
	if [ -d "${!1}" ]; then
		echo "exists."
	else
		echo "does not exist. Creating it."
		mkdir ${!1}
	fi
}


for filename in $sourcepath/*.dat; do

		[[ ${filename} =~ _T=([[:digit:]]*\.[[:digit:]]*)_ ]] && temperature=${BASH_REMATCH[1]}
		currentpath=${targetpath}/T=$temperature
		make_directory_if_necessary currentpath

		[[ ${filename} =~ _p=([[:digit:]]*\.[[:digit:]]*)_ ]] && pressure=${BASH_REMATCH[1]}
		currentpath=${currentpath}/p=$pressure
		make_directory_if_necessary currentpath

		mv $filename $currentpath

done
