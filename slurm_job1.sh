#!/bin/sh -f
#
###############################
# Specify your SLURM directives
###############################
#
# Job name:
#SBATCH --job-name BEM1
#SBATCH --output=slurm-%j-%x.out
#SBATCH --error=slurm-%j-%x.err
#
# If you want to specify a certain number of nodes
# and exactly 'ntasks-per-node' cpus on each node.
##SBATCH --nodes=4
##SBATCH --ntasks-per-node=24
#
# For hybrid mpi: e.g. 1 mpi process with
# n openmp threads
# SBATCH --ntasks=1
# SBATCH --cpus-per-task=1
# 
# Allocate full node and block for other jobs
#SBATCH --exclusive
#
# Walltime:
# SBATCH --time=72:00:00
###########################################

# Set OMP specific environment

module load openmpi/gcc/4.0.0

echo "Loaded mpi"
 
export OMP_NUM_THREADS=1

srun /home/bartsch/BEM/bem #omp_test/test_program
