#!/bin/bash

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=phi
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/data_analysis/avg_data_files"

sourcepath=$1
filename=$2
numberOfHeaderLines=$3
numberOfColumns=$4

submitFilename="computeAvgDataFilesFor_${filename}.sh"
srun_command="srun --ntasks=1 --error=./error_stream_${filename}.err ./compute_avg_data_files.sh ${sourcepath} "${filename}.dat" ${numberOfHeaderLines} ${numberOfColumns} &
wait"

echo "$settings" > ${submitFilename}
echo -e >> ${submitFilename}
echo "$srun_command" >> ${submitFilename}
sbatch ${submitFilename}
rm ${submitFilename}

