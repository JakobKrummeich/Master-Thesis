#!/bin/bash

settings="#!/bin/bash
#SBATCH --partition=normal,oip,phi
#SBATCH --nodes=1-1
#SBATCH --exclude=knightslanding02,knightslanding04
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/data_analysis/avg_data_files"

sourcepath=$1
filename=$2
numberOfHeaderLines=$3
numberOfColumns=$4

submitFilename="computeAvgDataFilesFor_${filename}.sh"
srun_command="srun ./compute_avg_data_files.sh ${sourcepath} "${filename}.dat" ${numberOfHeaderLines} ${numberOfColumns} &
wait"

echo "$settings" > ${submitFilename}
echo -e >> ${submitFilename}
echo "$srun_command" >> ${submitFilename}
sbatch ${submitFilename}
rm ${submitFilename}

