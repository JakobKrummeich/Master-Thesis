settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/semi_gcmc_no_field"

for N in {500,1000,2000,4000,8000,16000}; do

	for temperature in 0.7{500,495,490,485,480,475,470,465,460,455,450}; do
		echo "$settings" > N=${N}_T=${temperature}.sh
		echo  -e  >> N=${N}_T=${temperature}.sh
		srun_command="srun --ntasks=1 --error=error_stream_output/N=${N}_%J.err ./semi_gcmc_single_temp_${N} ${temperature} /home1/krummeich/Master-Thesis/data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh=0.75_deltaT=0.002/N=${N}/T=0.750000/States_N=${N}_T=0.750000_AvgDens=0.750000_MCRuns=300000_epsAB=0.100000.dat 10 &

wait"
		echo "$srun_command" >> N=${N}_T=${temperature}.sh
		sbatch N=${N}_T=${temperature}.sh
	done
done
