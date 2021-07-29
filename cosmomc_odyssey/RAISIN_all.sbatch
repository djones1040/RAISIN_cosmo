#!/bin/bash
#SBATCH --job-name=RAISIN_stat
#SBATCH --time=34:00:00
###SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=1-2
#SBATCH --cpus-per-task=1
#SBATCH --partition=broadwl-lc
#SBATCH --output=/scratch/midway2/rkessler/djones/cosmomc/chains/RAISIN_stat.log
#SBATCH --account=pi-rkessler
#SBATCH --mem=20GB

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module unload openmpi
module load intelmpi/5.1+intel-16.0
module load cfitsio/3
module load mkl
module load gcc/6.1
source /project2/rkessler/SURVEYS/PS1MD/USERS/djones/RAISIN/code/plc_3.0/plc-3.01/bin/clik_profile.sh

PARAMS=`expr ${SLURM_ARRAY_TASK_ID} - 1`

INI_FILES=(cosmomc/RAISIN_all.ini cosmomc/RAISIN_all.ini)
DONE_FILES=(done_0.txt done_1.txt)

cd /scratch/midway2/rkessler/djones/cosmomc/chains_2015/
mpirun /project2/rkessler/PRODUCTS/CosmoMC/v03/CosmoMC-master/cosmomc ${INI_FILES[$PARAMS]}
#mpirun $RAISIN_ROOT/CosmoMC/cosmomc ${INI_FILES[$PARAMS]}

if [ $? -eq 0 ]; then
    echo "SUCCESS" > ${DONE_FILES[$PARAMS]}
else
    echo "FAILURE" > ${DONE_FILES[$PARAMS]}
fi


