#!/bin/bash                                                                   
#SBATCH --nodes=1                                                             
#SBATCH --ntasks-per-node=4                                                   
#SBATCH --mem=230000MB                                                        
#SBATCH --gres=gpu:4                                                          
#SBATCH --time=00:30:00                                                       
#SBATCH --account=Sis21_baroni_0                                              
#SBATCH --partition=m100_usr_prod                                             
#SBATCH --job-name=NVT-$i                                                     
#SBATCH --mail-user=cmalosso@sissa.it                                         
#SBATCH --qos=m100_qos_dbg # da scommentare per la coda debug massimo 2 Nodi 
#SBATCH --mail-type=ALL                                                       


cd ${SLURM_SUBMIT_DIR}
echo $SLURM_SUBMIT_DIR
find  -type f  -exec touch {} +
