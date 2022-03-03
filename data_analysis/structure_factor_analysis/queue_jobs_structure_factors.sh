sourcepath=$1
targetFiles=$2
kAll=$3
kMax=$4
numberOfkIntervals=$5
numberOfAngleIntervals=$6
totalNumberOfParticles=$7

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=phi
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/data_analysis/structure_factor_analysis/"

filename="structureFactorComputation.sh"
echo "$settings" > $filename
echo  -e  >> $filename

srun_command="srun --ntasks=1 ./compute_structure_factors.sh ${sourcepath} ${targetFiles} ${kAll} ${kMax} ${numberOfkIntervals} ${numberOfAngleIntervals} ${totalNumberOfParticles} &

wait"

echo "$srun_command" >> $filename
sbatch $filename
rm $filename
