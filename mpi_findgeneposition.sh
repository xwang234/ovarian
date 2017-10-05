#!/usr/bin/env bash
#SBATCH -t 0-10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org
#SBATCH -n 10
hostname
echo "-----"
srun hostname
mpirun -n 1 mpi_findgeneposition.R
