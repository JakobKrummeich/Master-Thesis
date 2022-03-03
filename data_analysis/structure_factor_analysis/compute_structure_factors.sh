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
kAll=$3
kMax=$4
numberOfkIntervals=$5
numberOfAngleIntervals=$6
totalNumberOfParticles=$7

targetFilepath="${sourcepath}/aggregated_data/"
makeDirectoryIfNecessary targetFilepath

for runDir in ${sourcepath}/singleRunData/*; do

	files+=(${runDir}/${filename})
	
done

echo ${files[@]} | ./structure_factor "N=${totalNumberOfParticles}" ${targetFilepath} ${kAll} ${kMax} ${numberOfkIntervals} ${numberOfAngleIntervals} ${totalNumberOfParticles}

