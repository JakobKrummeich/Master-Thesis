settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/cmc_no_field"

declare -a equilibrium_values=(
	"0.750 0.5"
	"0.760 0.5"
	"0.770 0.5"
	"0.780 0.5"
	"0.790 0.5"
	"0.800 0.5"
	"0.810 0.5"
	"0.820 0.5"
	"0.830 0.5"
	"0.840 0.5"
)


for equilibrium_value in "${equilibrium_values[@]}"; do
	read -a tuple <<< "$equilibrium_value"
	filename="T=${tuple[0]}.sh"
	echo "$settings" > $filename
	echo  -e  >> $filename
	srun_command="srun --ntasks=1 --error=error_stream_output/T=${tuple[0]}_%J.err ./cmc ${tuple[0]} ${tuple[1]} 0 100000 0 20000000 &

wait"
	echo "$srun_command" >> $filename
	sbatch $filename
done

