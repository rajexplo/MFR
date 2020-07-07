#!/bin/bash
set -e
make clean
module load gcc/9.2.0
make -j8
sbatch HJB.sbatch 
make clean

