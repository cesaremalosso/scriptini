#!/bin/bash
for i in $1 
do
cat > NVTlammps0.in <<EOF
log "log0.log"
units metal                   
atom_style atomic             
read_data "$3"
pair_style deepmd $4
pair_coeff  

displace_atoms all random 0.02 0.02 0.02 4139$i units box
velocity all create $i 4154$i mom yes rot yes dist gaussian
                                                                   
timestep 0.2e-3
thermo 50                                                                                      

thermo_style custom step time ke pe etotal temp press density lx ly lz 

fix finvt1 all nvt temp $i $i \$(100.0*dt) 
dump NVTdump all custom 50 dump0.10fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass
dump binDump4 all custom 4 dump0.0.8fs.bin id type x y z vx vy vz
run $2 


write_restart final0.restart 
write_data final0.data
EOF


cat > slurm_simulation0.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=230000MB
#SBATCH --gpus-per-node=4
#SBATCH --time=05:00:00
#SBATCH --account=Sis21_baroni_0
#SBATCH --partition=m100_usr_prod
#SBATCH --job-name=NVT-$i
#SBATCH --mail-user=cmalosso@sissa.it
##SBATCH --qos=m100_qos_dbg # da scommentare per la coda debug massimo 2 Nodi
#SBATCH --mail-type=ALL

source /m100/home/userexternal/dtisi000/deepmd-1.3.3/moduli_lammps.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR

mpirun -np 4 /m100/home/userexternal/dtisi000/deepmd-1.3.3/lammps29Oct2020_dtisi/src/lmp_mpi < NVTlammps0.in > NVTlammps0.out
EOF

done
