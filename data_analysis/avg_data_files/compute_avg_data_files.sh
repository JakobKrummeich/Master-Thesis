#!/bin/bash

function makeDirectoryIfNecessary() {
	echo -n "${!1} "
	if [ -d "${!1}" ]; then
		echo "exists."
	else
		echo "does not exist. Creating it."
		mkdir ${!1}
	fi
}

sourcepath=$1
filename=$2
numberOfHeaderLines=$3
numberOfColumns=$4

targetFilepath="${sourcepath}/aggregated_data/"
makeDirectoryIfNecessary targetFilepath

for temperatureDir in ${sourcepath}/T=*; do

	for runDir in ${temperatureDir}/*; do

		files+=(${runDir}/${filename})
	
	done

	[[ ${temperatureDir} =~ /T=([[:digit:]]*.[[:digit:]]*) ]] && temperature=${BASH_REMATCH[1]}

	echo 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2
	echo "T=${temperature} analysis running." 1>&2
	echo "++++++++++++++++++++++++++++++++++++++++" 1>&2

	results=$(echo ${files[@]} | ./avg_data_files ${targetFilepath} ${filename} ${numberOfHeaderLines} ${numberOfColumns})

done
