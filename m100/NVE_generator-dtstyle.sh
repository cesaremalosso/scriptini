#!/bin/bash

function usage
{
   echo "         --first first run"
   echo "         --last final run"
   echo "         -t temperature of the system"
   echo "         -l number of step of the simulation"
   echo "         -m NN model"
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
        -t )              shift; temp=$1
                          ;;
        -l )              shift; length=$1
                          ;;
        -m )              shift; model=$1
                          ;;
        -h | --help )     usage
                          ;;
        * )               usage
    esac
    shift
done
if [ $run2 = $run1 ]
then for i  in $run1
do
j=$(awk -v k=$i 'BEGIN{print k-1'})
jj=$(awk -v k=$length 'BEGIN{print k-1'})
cat > NVElammps$i.in <<EOF
log "log$i.log"
units metal                   
atom_style atomic             
read_restart final$j.restart
pair_style deepmd $model
pair_coeff  

timestep 0.2e-3
thermo 100                                                                                               

compute Tp all temp
variable time equal step*dt
variable Temp equal c_Tp
variable epress equal c_thermo_press
variable s4 equal c_thermo_press[4]
variable s5 equal c_thermo_press[5]
variable s6 equal c_thermo_press[6]

thermo_style custom step time pe etotal lx press temp 

fix finve all nve     
run 1
unfix finve

fix finve2 all nve
dump myDumpforTrain all custom 4 dump$i.08fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass
dump_modify myDumpforTrain format float %20.15g
dump binDump4 all custom 4 dump$i.0.8fs.bin id type xu yu zu vx vy vz
print "&&&& START NVE"

fix 5 all print 2 "\${time} \${epress} \${Temp} \${s4} \${s5} \${s6}" file stress$i.0.4fs.out screen "no" title "time Press Temp sxy sxz syz" 

run $jj
unfix finve2
unfix 5

write_restart final$i.restart
write_data final$i.data
EOF

cat > slurm_simulation$i.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=230000MB
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --account=Sis21_baroni_0
#SBATCH --partition=m100_usr_prod
#SBATCH --job-name=NVE-$i
#SBATCH --mail-user=cmalosso@sissa.it
##SBATCH --qos=m100_qos_dbg # da scommentare per la coda debug massimo 2 Nodi
#SBATCH --mail-type=ALL

source /m100/home/userexternal/dtisi000/deepmd-1.3.3/moduli_lammps.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR

mpirun -np 4 /m100/home/userexternal/dtisi000/deepmd-1.3.3/lammps29Oct2020_dtisi/src/lmp_mpi < NVElammps$i.in > NVElammps$i.out  
EOF

done

else
for i in $(seq $run1 $run2)
do
j=$(awk -v k=$i 'BEGIN{print k-1'})
jj=$(awk -v k=$length 'BEGIN{print k-1'})
cat > NVElammps$i.in <<EOF
log "log$i.log"
units metal                   
atom_style atomic             
read_restart final$j.restart
pair_style deepmd $model
pair_coeff  

timestep 0.2e-3
thermo 100                                                                                               

compute Tp all temp
variable time equal step*dt
variable Temp equal c_Tp
variable epress equal c_thermo_press
variable s4 equal c_thermo_press[4]
variable s5 equal c_thermo_press[5]
variable s6 equal c_thermo_press[6]

thermo_style custom step time pe etotal lx press temp 

fix finve all nve     
run 1
unfix finve

fix finve2 all nve
dump myDumpforTrain all custom 4 dump$i.08fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass
dump_modify myDumpforTrain format float %20.15g
dump binDump4 all custom 4 dump$i.0.8fs.bin id type xu yu zu vx vy vz
print "&&&& START NVE"

fix 5 all print 2 "\${time} \${epress} \${Temp} \${s4} \${s5} \${s6}" file stress$i.0.4fs.out screen "no" title "time Press Temp sxy sxz syz" 

run $jj
unfix finve2
unfix 5
write_restart final$i.restart
write_data final$i.data
EOF
cat > slurm_simulation$i.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=230000MB
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --account=Sis21_baroni_0
#SBATCH --partition=m100_usr_prod
#SBATCH --job-name=NVE-$i
#SBATCH --mail-user=cmalosso@sissa.it
##SBATCH --qos=m100_qos_dbg # da scommentare per la coda debug massimo 2 Nodi
#SBATCH --mail-type=ALL

source /m100/home/userexternal/dtisi000/deepmd-1.3.3/moduli_lammps.sh

cd \${SLURM_SUBMIT_DIR}
echo \$SLURM_SUBMIT_DIR

mpirun -np 4 /m100/home/userexternal/dtisi000/deepmd-1.3.3/lammps29Oct2020_dtisi/src/lmp_mpi < NVElammps$i.in > NVElammps$i.out  
EOF
done
fi
