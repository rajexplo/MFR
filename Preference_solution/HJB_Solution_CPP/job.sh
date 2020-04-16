#!/bin/bash
set -e

make -j8
sbatch HJB.sbatch 
make clean

