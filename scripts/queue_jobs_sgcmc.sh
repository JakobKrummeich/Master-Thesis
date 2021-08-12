settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/semi_gcmc_no_field"

states_to_skip_per_run="1"
number_of_equilibration_sweeps="2000000"
max_number_of_sweeps="20000000"

for N in 16000; do
	for temperature in {0.7468,0.7469,0.7470,0.7471,0.7472}; do
		for run_offset in {0,20,40,60}; do
			submit_filename="N=${N}_T=${temperature}_${run_offset}.sh"
			echo "$settings" > ${submit_filename}
			echo  -e  >> ${submit_filename}
			srun_command="srun --ntasks=1 --error=error_stream_output/N=${N}_%J.err ./semi_gcmc_single_temp_${N} ${temperature} "/home1/krummeich/Master-Thesis/data/semi_gcmc_no_field/sorted_data/finite_size_scaling_roh=0.75_long_runs/N=${N}/T=0.750000/States_N=${N}_T=0.750000_AvgDens=0.750000_MCRuns=10000000_epsAB=0.100000.dat" ${states_to_skip_per_run} ${run_offset} ${number_of_equilibration_sweeps} ${max_number_of_sweeps} &

wait"
			echo "$srun_command" >> ${submit_filename}
			sbatch ${submit_filename}
		done
	done
done
