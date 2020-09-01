#!/bin/bash

#SBATCH --job-name=HJB_Solver
#SBATCH --error=job_%j.err
#SBATCH --output=job_%j.out
#SBATCH --partition=broadwl
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --qos=stafftest
#module load gcc/9.2.0

export OMP_NUM_THREADS=24

./HJBSolver ../TCRE_MacDougallEtAl2017_update.csv
