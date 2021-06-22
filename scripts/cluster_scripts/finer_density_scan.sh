#!/bin/bash 
#SBATCH --ntasks=8
#SBATCH --partition=test
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null
#SBATCH --mem-per-cpu=800mb

srun --ntasks=1 --error=output/AvgDens=0.51_%J.err ./semi_gcmc 0.51 0.66 0.016 0.35 1000000 &
srun --ntasks=1 --error=output/AvgDens=0.53_%J.err ./semi_gcmc 0.53 0.66 0.016 0.35 1000000 &
srun --ntasks=1 --error=output/AvgDens=0.55_%J.err ./semi_gcmc 0.55 0.66 0.016 0.35 1000000 &
srun --ntasks=1 --error=output/AvgDens=0.57_%J.err ./semi_gcmc 0.57 0.66 0.016 0.35 1000000 &
srun --ntasks=1 --error=output/AvgDens=0.59_%J.err ./semi_gcmc 0.59 0.66 0.016 0.35 1000000 &
srun --ntasks=1 --error=output/AvgDens=0.61_%J.err ./semi_gcmc 0.61 0.66 0.016 0.35 1000000 &
srun --ntasks=1 --error=output/AvgDens=0.63_%J.err ./semi_gcmc 0.63 0.66 0.016 0.35 1000000 &
srun --ntasks=1 --error=output/AvgDens=0.65_%J.err ./semi_gcmc 0.65 0.66 0.016 0.35 1000000 &


wait
