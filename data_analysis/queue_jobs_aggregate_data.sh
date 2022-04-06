#!/bin/bash

sourcepath=$1

cd avg_data_files
./analyze_singleRunData.sh ${sourcepath}

cd ../energy_analysis
./queue_jobs_energy_aggregation.sh ${sourcepath}

cd ../NA_analysis
./queue_jobs_NA_analysis.sh ${sourcepath}
