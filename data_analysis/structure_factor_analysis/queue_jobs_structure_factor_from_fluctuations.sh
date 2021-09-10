sourcepath=$1

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/data_analysis/structure_factor_analysis/"


submit_filename="structure_factors_from_fluctuations.sh"
echo "$settings" > ${submit_filename}

srun_command="srun --ntasks=1 --error=error_stream_output/structure_factors_from_fluctuations_%J.err ./compute_structure_factors_from_fluctuations.sh ${sourcepath} 0 0 &

wait"

echo "${srun_command}" >> ${submit_filename}
sbatch ${submit_filename}
