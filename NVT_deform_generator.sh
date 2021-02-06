#!/bin/bash
for i in $1 
do
cd $i'K'
cat > NVTlammps0.in <<EOF
log "log0.log"
units metal                   
atom_style atomic             
read_data "lammpsinit.data"
replicate 2 1 1 
pair_style deepmd frozen_model_SCAN_ultimo.pb
pair_coeff  

displace_atoms all random 0.02 0.02 0.02 4139$i units box
velocity all create $i 4154$i mom yes rot yes dist gaussian
                                                                   
timestep 0.2e-3
thermo 50                                                                                      

thermo_style custom step time ke pe etotal temp press density lx ly lz 


fix deform1 all deform 5 x final 0.0 19.549631828399999 y final 0.0 19.549631828399999 z  final 0.0 19.549631828399999  remap x 

fix finvt1 all nvt temp $i $i \$(100.0*dt) 
run 25000
unfix finvt1

fix finvt2 all nvt temp $i $i \$(100.0*dt) 
run 100000

dump NVTdump all custom 50 dump0.10fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass

write_restart final0.restart 
write_data final0.data
EOF


cat > slurm_simulation0.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=2
#SBATCH --time=11:59:00
#SBATCH --account=cmalosso
#SBATCH --partition=gpu1
#SBATCH --mem=63500
#SBATCH --job-name=$i.0NVT
#SBATCH --mail-user=cmalosso@sissa.it
#SBATCH --mail-type=ALL

source /home/dtisi/deepmd-tf2/moduli_anaconda.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=2
/opt/sissa/sissa-devtoolset7-openmpi-2.1.3/root/bin/mpirun -np 20 /home/dtisi/deepmd-tf2/lammps-master/src/lmp_mpi < NVTlammps0.in > NVTlammps0.out
EOF

cd ../
done
