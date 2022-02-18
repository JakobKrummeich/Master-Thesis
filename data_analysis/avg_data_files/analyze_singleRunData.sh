#!/bin/bash

sourcepath=$1

./queue_jobs_avg_data_files.sh ${sourcepath} "avgStresses_ShearStressSeries" "1" "1"
./queue_jobs_avg_data_files.sh ${sourcepath} "avgStresses_xEdges" "2" "5"
./queue_jobs_avg_data_files.sh ${sourcepath} "avgStresses_yEdges" "2" "5"
./queue_jobs_avg_data_files.sh ${sourcepath} "avgVelocities" "15" "0"
./queue_jobs_avg_data_files.sh ${sourcepath} "energySeries" "1" "3"
./queue_jobs_avg_data_files.sh ${sourcepath} "Img22" "1" "4"

