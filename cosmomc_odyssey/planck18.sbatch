#!/bin/bash
#SBATCH --job-name=RAISIN_stat
#SBATCH --time=34:00:00
#SBATCH --ntasks=8
#SBATCH --array=1-8
#SBATCH --cpus-per-task=1
#SBATCH --output=RAISIN_stat.log
#SBATCH --mem=30GB
#SBATCH -p shared

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

PARAMS=`expr ${SLURM_ARRAY_TASK_ID} - 1`

module load intel/2017c
source /n/holystore01/LABS/berger_lab/Lab/djones01/code/plc_3.0/plc-3.1/bin/clik_profile.sh

INI_FILES=(cosmomc_odyssey/planck18.ini cosmomc_odyssey/planck18.ini cosmomc_odyssey/planck18.ini cosmomc_odyssey/planck18.ini cosmomc_odyssey/planck18.ini cosmomc_odyssey/planck18.ini cosmomc_odyssey/planck18.ini cosmomc_odyssey/planck18.ini)
DONE_FILES=(done_0.txt done_2.txt done_3.txt done_4.txt done_5.txt done_6.txt done_7.txt done_8.txt)

cd /n/holystore01/LABS/berger_lab/Lab/djones01/RAISIN/chains
#mpirun /project2/rkessler/PRODUCTS/CosmoMC/v03/CosmoMC-master/cosmomc ${INI_FILES[$PARAMS]}
mpirun /n/holystore01/LABS/berger_lab/Lab/djones01/CosmoMC/cosmomc ${INI_FILES[$PARAMS]}

if [ $? -eq 0 ]; then
    echo "SUCCESS" > ${DONE_FILES[$PARAMS]}
else
    echo "FAILURE" > ${DONE_FILES[$PARAMS]}
fi
