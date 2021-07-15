settings="#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/cmc_no_field"

declare -a equilibrium_values=(
	"0.732 0.121744"
	"0.734 0.129572"
	"0.736 0.137313"
	"0.738 0.143692"
	"0.740 0.152967"
	"0.742 0.176043"
	"0.744 0.185197"
)


for equilibrium_value in "${equilibrium_values[@]}"; do
	read -a tuple <<< "$equilibrium_value"
	filename="T=${tuple[0]}.sh"
	echo "$settings" > $filename
	echo  -e  >> $filename
	srun_command="srun --ntasks=1 --error=error_stream_output/T=${tuple[0]}_%J.err ./cmc ${tuple[0]} ${tuple[1]} &

wait"
	echo "$srun_command" >> $filename
	sbatch $filename
done

