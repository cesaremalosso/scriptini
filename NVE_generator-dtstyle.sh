#!/bin/bash

function usage
{
   echo "         --first first run"
   echo "         --last final run"
   echo "         --gpu gpu for the job"
   echo "         -t temperature of the system"
   echo "         -l number of step of the simulation"
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
        -l )              shift; length=$1
                          ;;
        -h | --help )     usage
                          exit
                          ;;
        * )               usage
                          exit 
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
pair_style deepmd frozen_model.pb
pair_coeff  


timestep 0.2e-3

thermo_style custom step time pe etotal  press temp
                                      
                                                                   
thermo 100                                                                                               

compute TYPE all chunk/atom type
compute vcm all vcm/chunk TYPE                                                                                                         
compute energia_cinetica all ke/atom
compute energia_potenziale all pe/atom
compute flusso all heat/flux/full energia_cinetica energia_potenziale
compute Tp all temp

variable time equal step*dt
variable Temp equal c_Tp
variable flx1 equal c_flusso[1]
variable flx2 equal c_flusso[2]
variable flx3 equal c_flusso[3]
variable flc1 equal c_flusso[4]
variable flc2 equal c_flusso[5]
variable flc3 equal c_flusso[6]
variable flv1 equal c_flusso[7]
variable flv2 equal c_flusso[8]
variable flv3 equal c_flusso[9]
variable vcmOx equal c_vcm[1][1]
variable vcmOy equal c_vcm[1][2]
variable vcmOz equal c_vcm[1][3]
variable vcmHx equal c_vcm[2][1]
variable vcmHy equal c_vcm[2][2]
variable vcmHz equal c_vcm[2][3]  
variable epot equal c_thermo_pe
variable ekin equal c_kep
variable epress equal c_thermo_press
variable s4 equal c_thermo_press[4]
variable s5 equal c_thermo_press[5]
variable s6 equal c_thermo_press[6]


thermo_style custom step time ke pe etotal lx density press temp c_flusso[1] c_flusso[2] c_flusso[3] c_flusso[4] c_flusso[5] c_flusso[6] c_flusso[7] c_flusso[8] c_flusso[9]
#thermo_style custom step time pe etotal lx press temp 

fix finve all nve     

dump myDumpforTrain all custom 4 dump$i.08fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass
dump_modify myDumpforTrain format float %20.15g
dump binDump4 all custom 4 dump$i.0.8fs.bin id type xu yu zu vx vy vz
print "&&&& START NVE"

fix 4 all print 4 "\${time} \${epress} \${Temp} \${flx1} \${flx2} \${flx3} \${vcmOx} \${vcmOy} \${vcmOz} \${vcmHx} \${vcmHy} \${vcmHz}" file multi_currenti$i.0.8fs.out screen "no" title "time Press Temp c_flx[1] c_flx[2] c_flx[3] c_vcmO[1] c_vcmO[2] c_vcmO[3] c_vcmH[1] c_vcmH[2] c_vcmH[3]" 
fix 4444 all print 4 "\${Temp} \${flx1} \${flx2} \${flx3} \${flc1} \${flc2} \${flc3} \${flv1} \${flv2} \${flv3} \${vcmOx} \${vcmOy} \${vcmOz} " file multi_currenti_xcv_$i.0.8fs.out screen "no" title "Temp c_flx[1] c_flx[2] c_flx[3] c_flene[1] c_flene[2] c_flene[3] c_flvir[1] c_flvir[2] c_flvir[3] c_vcmO[1] c_vcmO[2] c_vcmO[3] " 
fix 5 all print 2 "\${time} \${epress} \${Temp} \${s4} \${s5} \${s6}" file stress$i.0.4fs.out screen "no" title "time Press Temp sxy sxz syz" 

run $length
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

else
for i in $(seq $run1 $run2)
do
j=$(awk -v k=$i 'BEGIN{print k-1'})
cat > NVElammps$i.in <<EOF
log "log$i.log"
units metal                   
atom_style atomic             
read_restart final$j.restart
pair_style deepmd frozen_model.pb
pair_coeff  


timestep 0.2e-3

thermo_style custom step time pe etotal  press temp
                                      
                                                                   
thermo 100                                                                                               

compute TYPE all chunk/atom type
compute vcm all vcm/chunk TYPE                                                                                                         
compute energia_cinetica all ke/atom
compute energia_potenziale all pe/atom
compute flusso all heat/flux/full energia_cinetica energia_potenziale
compute Tp all temp

variable time equal step*dt
variable Temp equal c_Tp
variable flx1 equal c_flusso[1]
variable flx2 equal c_flusso[2]
variable flx3 equal c_flusso[3]
variable flc1 equal c_flusso[4]
variable flc2 equal c_flusso[5]
variable flc3 equal c_flusso[6]
variable flv1 equal c_flusso[7]
variable flv2 equal c_flusso[8]
variable flv3 equal c_flusso[9]
variable vcmOx equal c_vcm[1][1]
variable vcmOy equal c_vcm[1][2]
variable vcmOz equal c_vcm[1][3]
variable vcmHx equal c_vcm[2][1]
variable vcmHy equal c_vcm[2][2]
variable vcmHz equal c_vcm[2][3]  
variable epot equal c_thermo_pe
variable ekin equal c_kep
variable epress equal c_thermo_press
variable s4 equal c_thermo_press[4]
variable s5 equal c_thermo_press[5]
variable s6 equal c_thermo_press[6]


thermo_style custom step time ke pe etotal lx density press temp c_flusso[1] c_flusso[2] c_flusso[3] c_flusso[4] c_flusso[5] c_flusso[6] c_flusso[7] c_flusso[8] c_flusso[9]
#thermo_style custom step time pe etotal lx press temp 

fix finve all nve     

dump myDumpforTrain all custom 4 dump$i.08fs.lammpstrj id type x y z ix iy iz vx vy vz fx fy fz mass
dump_modify myDumpforTrain format float %20.15g
dump binDump4 all custom 4 dump$i.0.8fs.bin id type xu yu zu vx vy vz
print "&&&& START NVE"

fix 4 all print 4 "\${time} \${epress} \${Temp} \${flx1} \${flx2} \${flx3} \${vcmOx} \${vcmOy} \${vcmOz} \${vcmHx} \${vcmHy} \${vcmHz}" file multi_currenti$i.0.8fs.out screen "no" title "time Press Temp c_flx[1] c_flx[2] c_flx[3] c_vcmO[1] c_vcmO[2] c_vcmO[3] c_vcmH[1] c_vcmH[2] c_vcmH[3]" 
fix 4444 all print 4 "\${Temp} \${flx1} \${flx2} \${flx3} \${flc1} \${flc2} \${flc3} \${flv1} \${flv2} \${flv3} \${vcmOx} \${vcmOy} \${vcmOz} " file multi_currenti_xcv_$i.0.8fs.out screen "no" title "Temp c_flx[1] c_flx[2] c_flx[3] c_flene[1] c_flene[2] c_flene[3] c_flvir[1] c_flvir[2] c_flvir[3] c_vcmO[1] c_vcmO[2] c_vcmO[3] " 
fix 5 all print 2 "\${time} \${epress} \${Temp} \${s4} \${s5} \${s6}" file stress$i.0.4fs.out screen "no" title "time Press Temp sxy sxz syz" 

run $length
write_restart final$i.restart
write_data final$i.data
EOF
if [ "$gpu" = "gpu1" ]
then cat > slurm_simulation$i.pbs <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
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
export OMP_NUM_THREADS=20
/opt/sissa/sissa-devtoolset7-openmpi-2.1.3/root/bin/mpirun -np 1 /home/dtisi/deepmd-tf2/lammps-master/src/lmp_mpi < NVElammps$i.in > NVElammps$i.out
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
