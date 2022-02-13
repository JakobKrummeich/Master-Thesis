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

for runDir in ${sourcepath}/singleRunData/*; do

	files+=(${runDir}/${filename})
	
done


results=$(echo ${files[@]} | ./avg_data_files ${targetFilepath} ${filename} ${numberOfHeaderLines} ${numberOfColumns})


