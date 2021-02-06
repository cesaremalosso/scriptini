#!/bin/bash

function usage
{
   echo "         --first first run"
   echo "         --last final run"
   echo "         --gpu gpu for the job"
   echo "         -t temperature of the system"
   echo "         -h  --help      print this help"
   echo "Example:  NVE_generator.sh -t 500 --first 1 --last 10"
}
run2=0
while [ $# -gt 0 ]; do
    case $1 in
        --first )         shift; run1=$1
                          ;;
        --last )         shift; run2=$1
                          ;;
        --gpu )           shift; gpu=$1
                          ;;
        -t )              shift; temp=$1
                          ;;
        -h | --help )     usage
                          exit
                          ;;
        * )               usage
                          exit 1
    esac
    shift
done
if [ $run2 = $run1 ]
then for i  in $run1
do
j=$(awk -v k=$i 'BEGIN{print k-1'})
cat > NVElammps$i.in <<EOF
log "log$i.log"
units metal                   
atom_style atomic             
read_restart final$j.restart
pair_style deepmd frozen_model_SCAN_ultimo.pb
pair_coeff  

timestep 0.2e-3

                                      
thermo 50            
                                                                                   
variable t equal time
variable epot equal pe
variable etot equal etotal
variable temp equal temp


variable stressxx equal c_thermo_press[1]
variable stressyy equal c_thermo_press[2]
variable stresszz equal c_thermo_press[3]

variable stressxy equal c_thermo_press[4]
variable stressxz equal c_thermo_press[5]
variable stressyz equal c_thermo_press[6]

thermo_style custom step time ke pe etotal temp press density 
fix finve all nve     

dump mydump all custom 50 dump$i.10fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass
dump bindump all custom 5 dump$i.1fs.bin id type xu yu zu vx vy vz
print "&&&& START NVE"

fix myfix all print 2 '\${t} \${etot} \${epot} \${temp} \${stressxy} \${stressxz} \${stressyz}' file offstress$i.0.4fs.out screen "no" title "time etot epot temp stress_xy stress_xz stress_yz"
fix myfix1 all print 2 '\${t} \${etot} \${epot} \${temp} \${stressxx} \${stressyy} \${stresszz}' file instress$i.0.4fs.out screen "no" title "time etot epot temp stress_xx stress_yy stress_zz"

run 125000 
write_restart final$i.restart
write_data final$i.data
EOF

if [ "$gpu" = "gpu1" ]
then cat > slurm_simulation$i.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=2
#SBATCH --time=11:59:00
#SBATCH --account=cmalosso
#SBATCH --partition=$gpu
#SBATCH --mem=63500
#SBATCH --job-name=$temp.0NVE
#SBATCH --mail-user=cmalosso@sissa.it
#SBATCH --mail-type=ALL

source /home/dtisi/deepmd-tf2/moduli_anaconda.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=2
/opt/sissa/sissa-devtoolset7-openmpi-2.1.3/root/bin/mpirun -np 20 /home/dtisi/deepmd-tf2/lammps-master/src/lmp_mpi < NVElammps$i.in > NVElammps$i.out
EOF

elif [ "$gpu" = "gpu2" ]
then cat > slurm_simulation$i.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=63500
#SBATCH --time=11:59:00
#SBATCH --account=cmalosso
#SBATCH --partition=$gpu
#SBATCH --job-name=$temp.0NVE
#SBATCH --mail-user=cmalosso@sissa.it
#SBATCH --mail-type=ALL

source /home/dtisi/deepmd-tf2/moduli_anaconda.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=2
/opt/sissa/sissa-devtoolset7-openmpi-2.1.3/root/bin/mpirun -np 1 /home/dtisi/deepmd-tf2/lammps-master/src/lmp_mpi < NVElammps$i.in > NVElammps$i.out
EOF
else
echo '--gpus value not admitted'
fi
done

else
for i in $(seq $run1 $run2)
do
j=$(awk -v k=$i 'BEGIN{print k-1'})
cat > NVElammps$i.in <<EOF
log "log$i.log"
units metal                   
atom_style atomic             
read_restart final$j.restart
pair_style deepmd frozen_model_SCAN_ultimo.pb
pair_coeff  

timestep 0.2e-3

                                      
thermo 50            
                                                                                   
variable t equal time
variable epot equal pe
variable etot equal etotal
variable temp equal temp


variable stressxx equal c_thermo_press[1]
variable stressyy equal c_thermo_press[2]
variable stresszz equal c_thermo_press[3]

variable stressxy equal c_thermo_press[4]
variable stressxz equal c_thermo_press[5]
variable stressyz equal c_thermo_press[6]

thermo_style custom step time ke pe etotal temp press density 
fix finve all nve     

dump mydump all custom 50 dump$i.10fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass
dump bindump all custom 5 dump$i.1fs.bin id type xu yu zu vx vy vz
print "&&&& START NVE"

fix myfix all print 2 '\${t} \${etot} \${epot} \${temp} \${stressxy} \${stressxz} \${stressyz}' file offstress$i.0.4fs.out screen "no" title "time etot epot temp stress_xy stress_xz stress_yz"
fix myfix1 all print 2 '\${t} \${etot} \${epot} \${temp} \${stressxx} \${stressyy} \${stresszz}' file instress$i.0.4fs.out screen "no" title "time etot epot temp stress_xx stress_yy stress_zz"

run 125000 
write_restart final$i.restart
write_data final$i.data
EOF
if [ "$gpu" = "gpu1" ]
then cat > slurm_simulation$i.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=2
#SBATCH --time=11:59:00
#SBATCH --account=cmalosso
#SBATCH --partition=$gpu
#SBATCH --mem=63500
#SBATCH --job-name=$temp.0NVE
#SBATCH --mail-user=cmalosso@sissa.it
#SBATCH --mail-type=ALL

source /home/dtisi/deepmd-tf2/moduli_anaconda.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=2
/opt/sissa/sissa-devtoolset7-openmpi-2.1.3/root/bin/mpirun -np 20 /home/dtisi/deepmd-tf2/lammps-master/src/lmp_mpi < NVElammps$i.in > NVElammps$i.out
EOF

elif [ "$gpu" = "gpu2" ]
then cat > slurm_simulation$i.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=63500
#SBATCH --time=11:59:00
#SBATCH --account=cmalosso
#SBATCH --partition=$gpu
#SBATCH --job-name=$temp.0NVE
#SBATCH --mail-user=cmalosso@sissa.it
#SBATCH --mail-type=ALL

source /home/dtisi/deepmd-tf2/moduli_anaconda.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=36
/opt/sissa/sissa-devtoolset7-openmpi-2.1.3/root/bin/mpirun -np 1 /home/dtisi/deepmd-tf2/lammps-master/src/lmp_mpi < NVElammps$i.in > NVElammps$i.out
EOF
else
echo '--gpus value not admitted'
fi
done
fi
