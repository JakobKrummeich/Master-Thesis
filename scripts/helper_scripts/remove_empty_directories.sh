#!/bin/bash

targetDir=$1

for dir in ${targetDir}*/ ; do

	if [ -d "${dir}" ]; then
		if [ ! "$(ls -A ${dir})" ]; then

			rm -r ${dir}
			echo "Removing ${dir}."

		fi
	fi

done


