#!/bin/bash

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=phi
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/data_analysis/NA_analysis"

sourcepath=$1

submitFilename="computeNA_analysis.sh"
srun_command="srun --ntasks=1 --error=./error_stream.err ./compute_aggregated_distribution.sh ${sourcepath} &
wait"

echo "$settings" > ${submitFilename}
echo -e >> ${submitFilename}
echo "$srun_command" >> ${submitFilename}
sbatch ${submitFilename}
rm ${submitFilename}

