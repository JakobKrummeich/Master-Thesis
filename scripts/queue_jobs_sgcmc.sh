settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/semi_gcmc_no_field"

for N in {500,1000,2000,4000,8000,16000}; do

	for temperature in {0.74680,0.74690,0.74700,0.74710,0.74720}; do
		for run_offset in {0,20,40,60}; do
			echo "$settings" > N=${N}_T=${temperature}_${run_offset}.sh
			echo  -e  >> N=${N}_T=${temperature}_${run_offset}.sh
			srun_command="srun --ntasks=1 --error=error_stream_output/N=${N}_%J.err ./semi_gcmc_single_temp_${N} ${temperature} /home1/krummeich/Master-Thesis/data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh=0.75_long_runs/N=${N}/T=0.750000/States_N=${N}_T=0.750000_AvgDens=0.750000_MCRuns=10000000_epsAB=0.100000.dat 10 ${run_offset} &

wait"
			echo "$srun_command" >> N=${N}_T=${temperature}_${run_offset}.sh
			sbatch N=${N}_T=${temperature}_${run_offset}.sh
		done
	done
done
