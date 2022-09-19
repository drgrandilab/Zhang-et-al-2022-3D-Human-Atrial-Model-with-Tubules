#!/bin/bash
#SBATCH --job-name=T10_1hz
#SBATCH --output=res_omp.txt
#SBATCH --error=err_omp.txt
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

make clean; make; time ./pace_2 1000 28 "pool_tubule/tub_input_ver2_10.txt"