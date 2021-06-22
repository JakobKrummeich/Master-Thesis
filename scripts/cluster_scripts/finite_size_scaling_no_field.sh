#!/bin/bash 
#SBATCH --ntasks=1
#SBATCH --partition=parallel
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/simulation/semi_gcmc_no_field

srun --ntasks=1 --error=error_stream_output/N=16000_%J.err ./semi_gcmc_16000 0.6 0.62 0.002 0.58 300000 0.1 20 20 ../../data/semi_gcmc_no_field/fresh_data/N=16000 &
srun --ntasks=1 --error=error_stream_output/N=8000_%J.err  ./semi_gcmc_8000  0.6 0.62 0.002 0.58 300000 0.1 20 20 ../../data/semi_gcmc_no_field/fresh_data/N=8000 &
srun --ntasks=1 --error=error_stream_output/N=4000_%J.err  ./semi_gcmc_4000  0.6 0.62 0.002 0.58 300000 0.1 20 20 ../../data/semi_gcmc_no_field/fresh_data/N=4000 &
srun --ntasks=1 --error=error_stream_output/N=2000_%J.err  ./semi_gcmc_2000  0.6 0.62 0.002 0.58 300000 0.1 20 20 ../../data/semi_gcmc_no_field/fresh_data/N=2000 &
srun --ntasks=1 --error=error_stream_output/N=1000_%J.err  ./semi_gcmc_1000  0.6 0.62 0.002 0.58 300000 0.1 20 20 ../../data/semi_gcmc_no_field/fresh_data/N=1000 &
srun --ntasks=1 --error=error_stream_output/N=500_%J.err   ./semi_gcmc_500   0.6 0.62 0.002 0.58 300000 0.1 20 20 ../../data/semi_gcmc_no_field/fresh_data/N=500 &

wait
