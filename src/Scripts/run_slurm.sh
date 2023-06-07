#!/bin/bash

#SBATCH -n 8 -N 1 --exclusive -w node[014] --threads-per-core=1 -p mixedp  --error=slurm-%j.err --output=slurm-%j.out

export OMP_NUM_THREADS=1
ulimit -s unlimited

mpirun -n 8  ../exec/main
