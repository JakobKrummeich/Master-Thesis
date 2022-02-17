#!/bin/bash

sourcepath=$1

./queue_jobs_avg_data_files.sh ${sourcepath} "avgStresses_ShearStressSeries.dat" "1" "1"
./queue_jobs_avg_data_files.sh ${sourcepath} "avgStresses_xEdges.dat" "2" "5"
./queue_jobs_avg_data_files.sh ${sourcepath} "avgVelocities.dat" "15" "0"
./queue_jobs_avg_data_files.sh ${sourcepath} "energySeries.dat" "1" "3"
./queue_jobs_avg_data_files.sh ${sourcepath} "Img22.dat" "1" "4"

