#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=182000
#SBATCH --time=01:59:00
#SBATCH --account=Sis22_baroni_1
#SBATCH --partition=skl_usr_prod
#SBATCH --job-name=ttNVT
#SBATCH --mail-user=cmalosso@sissa.it
#SBATCH --mail-type=ALL


cd ${SLURM_SUBMIT_DIR}
echo $SLURM_SUBMIT_DIR
find  -type f  -exec touch {} +
