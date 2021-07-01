#!/bin/bash

sourcepath=$1
targetpath=$2
lower_limit=$3
upper_limit=$4

files=(${sourcepath}/*.dat)

for file in ${files[@]}; do

	[[ ${file} =~ epsAB=0\.100000_([[:digit:]]*) ]] && runnumber=${BASH_REMATCH[1]}

	echo $runnumber

	if (( runnumber > $lower_limit && runnumber <= $upper_limit )); then
		mv $file $targetpath
	fi
done
