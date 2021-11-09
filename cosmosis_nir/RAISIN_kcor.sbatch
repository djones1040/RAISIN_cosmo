#!/bin/bash -l
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --time=34:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=5000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --job-name=RAISIN_stat

mpirun -n 2 cosmosis --mpi RAISIN_kcor.ini

