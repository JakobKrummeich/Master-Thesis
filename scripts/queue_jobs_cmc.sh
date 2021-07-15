settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/cmc_no_field"

for temperature in  {0.84,0.83,0.82,0.81,0.80,0.79,0.78,0.77,0.76,0.75}; do
	echo "$settings" > T=${temperature}.sh
	echo  -e  >> T=${temperature}.sh
	srun_command="srun --ntasks=1 --error=error_stream_output/T=${temperature}_%J.err ./cmc ${temperature} 0.5 &

wait"
	echo "$srun_command" >> T=${temperature}.sh
	sbatch T=${temperature}.sh
done

