#!/bin/bash

sourcepath=$1

for filename in ${sourcepath}*/NA_Series*; do

		[[ ${filename} =~ _N=([[:digit:]]*)_ ]] && number_of_particles=${BASH_REMATCH[1]}
		[[ ${filename} =~ _T=([[:digit:]]*.[[:digit:]]*)_ ]] && temperature=${BASH_REMATCH[1]}
		
		./compute_first_moments $filename $number_of_particles

done




