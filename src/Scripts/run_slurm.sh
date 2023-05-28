#!/bin/bash

#SBATCH -n $2 -N $3 --exclusive -x node[015,024,041,006,013,028,029,007-012,014,016-023,026-027,038,043-046] --threads-per-core=1 -p mixedp  --error=slurm-%j.err --output=slurm-%j.out

## 30 - node[015,024,007,009,011,013,025,028,029,006,008,010,014,016-023,026-027,038,043-046]
## mantis - node[015,024,012,041,007,009,011,013,025,028,029,006,008,010,014,026-027,030-039,043-046]

export OMP_NUM_THREADS=1
ulimit -s unlimited
# module load intel/18.0.1.163 openmpi/4.1.3-intel

#load_oneapi
#module load intel/oneAPI/compiler/latest && source /opt/sw/intel-oneApi/hpc/setvars.sh

<LOAD_MODULES>
#module load intel/oneAPI/vtune/latest

#mpirun -n $2 vtune -collect hotspots -r res$2  ../exec/main
mpirun -n $2  ../exec/main
#mpirun -n 64 ./a.out
