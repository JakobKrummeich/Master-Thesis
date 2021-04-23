#!/bin/bash 
#SBATCH --ntasks=6
#SBATCH --partition=test
#SBATCH --nodes=1-1
#SBATCH --output=/dev/null

srun --ntasks=1 --error=output/AvgDens=0.80_%J.err ./semi_gcmc 0.80 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.75_%J.err ./semi_gcmc 0.75 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.70_%J.err ./semi_gcmc 0.70 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.65_%J.err ./semi_gcmc 0.65 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.60_%J.err ./semi_gcmc 0.60 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.55_%J.err ./semi_gcmc 0.55 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.50_%J.err ./semi_gcmc 0.50 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.45_%J.err ./semi_gcmc 0.45 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.40_%J.err ./semi_gcmc 0.40 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.35_%J.err ./semi_gcmc 0.35 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.30_%J.err ./semi_gcmc 0.30 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.25_%J.err ./semi_gcmc 0.25 0.9 0.25 0.8 10000 &
srun --ntasks=1 --error=output/AvgDens=0.20_%J.err ./semi_gcmc 0.20 0.9 0.25 0.8 10000 &

wait
