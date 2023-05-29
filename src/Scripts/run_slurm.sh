#!/bin/bash

#SBATCH -n 1 -N 1 --exclusive -x node[015,024,041,006,013,028,029,007-012,014,016-023,026-027,038,043-046] --threads-per-core=1 -p mixedp  --error=slurm-%j.err --output=slurm-%j.out

export OMP_NUM_THREADS=1
ulimit -s unlimited

mpirun -n 1  ../exec/main
