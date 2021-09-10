settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/semi_gcmc_no_field"

states_to_skip_per_run="1"
number_of_equilibration_sweeps="0"
max_number_of_sweeps="20000000"

for N in {16000,8000,4000,2000,1000,500}; do
	for temperature in {0.746800,0.746900,0.747000,0.747100,0.747200}; do
		for run_offset in {0,20,40,60}; do
			submit_filename="N=${N}_T=${temperature}_${run_offset}.sh"
			echo "$settings" > ${submit_filename}
			echo  -e  >> ${submit_filename}
			srun_command="srun --ntasks=1 --error=error_stream_output/N=${N}_%J.err ./semi_gcmc_single_temp_${N} ${temperature} "/home1/krummeich/Master-Thesis/data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh=0.75_long_runs_2/N=${N}/T=${temperature}/States_N=${N}_T=${temperature}_AvgDens=0.750000_MCRuns=20000000_epsAB=0.100000.dat" ${states_to_skip_per_run} ${run_offset} ${run_offset} ${number_of_equilibration_sweeps} ${max_number_of_sweeps} &

wait"
			echo "$srun_command" >> ${submit_filename}
			sbatch ${submit_filename}
		done
	done
done
