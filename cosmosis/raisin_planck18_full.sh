#!/bin/bash -l
###SBATCH -n 50 # Number of cores requested
###SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --time=34:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --ntasks=50
#SBATCH --array=1-50
#SBATCH --cpus-per-task=1
#SBATCH --mem=20000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname
#SBATCH --job-name=RAISIN_stat

export VALUESINCLUDE=blank_include.ini
export PRIORSINCLUDE=blank_include.ini

mpirun -n 50 cosmosis --mpi raisin_planck_full.ini
