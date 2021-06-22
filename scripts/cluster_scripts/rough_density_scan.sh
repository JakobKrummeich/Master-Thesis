#!/bin/bash 
#SBATCH --ntasks=8
#SBATCH --partition=test
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=800mb

srun --ntasks=1 --error=output/AvgDens=0.80_%J.err ./semi_gcmc 0.80 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.76_%J.err ./semi_gcmc 0.76 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.72_%J.err ./semi_gcmc 0.72 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.68_%J.err ./semi_gcmc 0.68 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.64_%J.err ./semi_gcmc 0.64 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.60_%J.err ./semi_gcmc 0.60 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.56_%J.err ./semi_gcmc 0.56 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.52_%J.err ./semi_gcmc 0.52 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.48_%J.err ./semi_gcmc 0.48 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.44_%J.err ./semi_gcmc 0.44 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.40_%J.err ./semi_gcmc 0.40 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.36_%J.err ./semi_gcmc 0.36 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.32_%J.err ./semi_gcmc 0.32 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.28_%J.err ./semi_gcmc 0.28 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.24_%J.err ./semi_gcmc 0.24 0.9 0.05 0.25 600000 &
srun --ntasks=1 --error=output/AvgDens=0.20_%J.err ./semi_gcmc 0.20 0.9 0.05 0.25 600000 &


wait
