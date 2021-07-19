#!/bin/bash
#SBATCH --job-name=RAISIN_cosmosis
#SBATCH --time=34:00:00
###SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=1-8
#SBATCH --cpus-per-task=1
#SBATCH --partition=broadwl-lc
#SBATCH --output=/scratch/midway2/rkessler/djones/cosmomc/chains/RAISIN_cosmosis_planck18.log
#SBATCH --account=pi-rkessler
#SBATCH --mem=30GB

export RUNNAME=des_sne_wcdm
export DATAFILE=des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits
export DATAFILE_SR=des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate_sr.npy
export SCALEFILE=blanck_include.ini
export INCLUDEFILE=fiducial/params_raisin_sne.ini
export VALUESINCLUDE=fiducial/values_w.ini
export PRIORSINCLUDE=blank_include.ini

mpirun -n 8 cosmosis --mpi fiducial/params.ini -p runtime.sampler='importance'
#cosmosis fiducial/params.ini
