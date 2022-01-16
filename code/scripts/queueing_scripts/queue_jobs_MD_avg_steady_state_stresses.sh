function make_directory_if_necessary() {
	echo -n "${!1} "
	if [ -d "${!1}" ]; then
		echo "exists."
	else 
		echo "does not exist. Creating it."
		mkdir ${!1}
	fi
}

temperature=2.0

datapath="data/MDSimulations/avg_steady_state_stresses/T=${temperature}" # relative to root directory of repository

new_directory="../../../${datapath}"

make_directory_if_necessary new_directory

settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/code/MDSimulation/avg_steady_state_stresses"

declare -a shear_rate_eq_steps_pairs=(
	"0.01 200000"
	"0.02 100000"
	"0.03 66667"
	"0.04 50000"
	"0.05 40000"
	"0.06 33333"
	"0.07 28571"
	"0.08 25000"
	"0.09 22222"
	"0.10 20000"
)

for shear_rate_eq_steps_pair in "${shear_rate_eq_steps_pairs[@]}"; do
	read -a tuple <<< "$shear_rate_eq_steps_pair"
	submit_filename="avg_stresses_N=1000_T=2.0_shear_rate=${tuple[0]}.sh"
	echo "$settings" > ${submit_filename}
	echo  -e  >> ${submit_filename}

	shear_directory=${new_directory}/shear_rate=${tuple[0]}
	make_directory_if_necessary shear_directory

	srun_command="srun --ntasks=1 --error=error_stream_output/N=${N}_%J.err ./avg_steady_state_stresses ${tuple[0]} ${tuple[1]} ${shear_directory} &

wait"
	echo "$srun_command" >> ${submit_filename}
	#sbatch ${submit_filename}
done

