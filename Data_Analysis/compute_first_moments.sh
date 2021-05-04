#!/bin/bash

sourcepath=$1

first_line_written=false

for filename in ${sourcepath}*/NA_Series*; do

		[[ ${filename} =~ _N=([[:digit:]]*)_ ]] && number_of_particles=${BASH_REMATCH[1]}
		[[ ${filename} =~ _T=([[:digit:]]*.[[:digit:]]*)_ ]] && temperature=${BASH_REMATCH[1]}
		[[ ${filename} =~ ('N='[[:digit:]]*).*_('AvgDens'.*) ]] && output_filename=phase_diagram_${BASH_REMATCH[1]}_${BASH_REMATCH[2]}

		if [ "$first_line_written" = false ] ; then
			echo -e  "temperature\txA" >> analysis_results/${output_filename}
			first_line_written=true;
		fi

		echo -n -e "${temperature}\t" >> analysis_results/${output_filename}
		./compute_first_moments $filename $number_of_particles >> analysis_results/${output_filename}

done




