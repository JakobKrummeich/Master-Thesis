#!/bin/bash

settings="#!/bin/bash
#SBATCH --partition=normal,oip,phi
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --exclude=knightslanding02,knightslanding04
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/data_analysis/NA_analysis"

sourcepath=$1

submitFilename="computeNA_analysis.sh"
srun_command="srun --ntasks=1 --ntasks-per-core=1 --error=./error_stream.err ./compute_aggregated_distribution.sh ${sourcepath} &
wait"

echo "$settings" > ${submitFilename}
echo -e >> ${submitFilename}
echo "$srun_command" >> ${submitFilename}
sbatch ${submitFilename}
rm ${submitFilename}

