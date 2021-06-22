#!/bin/bash 
#SBATCH --ntasks=8
#SBATCH --partition=test
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=800mb
#SBATCH --workdir=/home1/krummeich/simulation

srun --ntasks=1 --error=error_stream_output/eps=0.45_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.45 &
srun --ntasks=1 --error=error_stream_output/eps=0.40_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.40 &
srun --ntasks=1 --error=error_stream_output/eps=0.35_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.35 &
srun --ntasks=1 --error=error_stream_output/eps=0.30_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.30 &
srun --ntasks=1 --error=error_stream_output/eps=0.25_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.25 &
srun --ntasks=1 --error=error_stream_output/eps=0.20_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.20 &
srun --ntasks=1 --error=error_stream_output/eps=0.15_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.15 &
srun --ntasks=1 --error=error_stream_output/eps=0.10_%J.err ./semi_gcmc 0.6 1.0 0.035 0.31 500000 0.10 &


wait
