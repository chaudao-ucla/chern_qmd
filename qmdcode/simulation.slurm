#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cnd28@cornell.edu
#SBATCH --partition=standard
module load gcc
module load armadillo
./timd