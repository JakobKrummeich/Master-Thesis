#!/bin/bash 
#SBATCH --ntasks=20
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/Master-Thesis/simulation/pressure_mc

srun --ntasks=1 --error=error_stream_output/T=1.000_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 1.000 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.975_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.975 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.950_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.950 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.925_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.925 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.900_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.900 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.875_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.875 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.850_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.850 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.825_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.825 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.800_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.800 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.775_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.775 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.750_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.750 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.725_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.725 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.700_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.700 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.675_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.675 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.650_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.650 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.625_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.625 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.600_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.600 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.575_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.575 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.550_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.550 ../../data/pressure_mc/fresh_data &
srun --ntasks=1 --error=error_stream_output/T=0.525_%J.err ./pressure_mc 0.93 4.0 0.05 1.0 30000 0.1 1 0 0.525 ../../data/pressure_mc/fresh_data &

wait

