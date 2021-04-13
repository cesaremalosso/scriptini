import sys
sys.path.append("/home/cmalosso/scriptini")
import read_lammps_dump as rd
a_B = 0.529177210903
if len(sys.argv)<6:
    print("Usage: {} <dump_name_input> <POSCAR_name_output_directory> <every_n> <type_Li> <type_Cl> <type_O>".format(sys.argv[0]))
else:
    name_input=sys.argv[1]
    print(name_input)
    name_output=sys.argv[2]
    every_n=int(sys.argv[3])
    Li_type=int(sys.argv[4])
    Cl_type=int(sys.argv[5])
    O_type=int(sys.argv[6])
    traj = rd.LAMMPS_Dump(name_input)
    number_of_POSCARS=traj.TOT_TIMESTEPS//every_n
    traj.read_timesteps((traj.FIRST_TIMESTEP,traj.LAST_TIMESTEP+traj.DELTA_TIMESTEP,traj.DELTA_TIMESTEP*every_n))
    for nn in range(number_of_POSCARS):
        with open("{}/POSCAR.{}".format(name_output, nn), 'w') as out:
            out.write('POSCAR generated with dump2init \n')
            out.write('1.0 \n')
            #if (units == 'bohr'):
            #    out.write(str(a_B)+'\n')
            #else:
            #    print('only bohr to metal unit implemented')
            out.write('{}   {:.9f}   {:.9f}\n'.format(traj.BOX_BOUNDS[0,1],0,0))
            out.write('{:.9f}   {}   {:.9f}\n'.format(0,traj.BOX_BOUNDS[0,1],0))
            out.write('{:.9f}   {:.9f}   {}\n'.format(0,0,traj.BOX_BOUNDS[0,1]))
            out.write('Li      Cl      O\n')
            out.write('{}   {}   {}\n'.format(80, 26, 27))
            out.write('Cartesian\n')
            for n in range(traj.NATOMS):
                if (traj.data[nn]['type'][n]==Li_type):
                    out.write('{}   {}   {}\n'.format(traj.data[nn]['xu'][n,0],traj.data[nn]['yu'][n,0],traj.data[nn]['zu'][n,0]))
            for n in range(traj.NATOMS):
                if (traj.data[nn]['type'][n]==Cl_type):
                    out.write('{}   {}   {}\n'.format(traj.data[nn]['xu'][n,0],traj.data[nn]['yu'][n,0],traj.data[nn]['zu'][n,0]))
            for n in range(traj.NATOMS):
                if (traj.data[nn]['type'][n]==O_type):
                    out.write('{}   {}   {}\n'.format(traj.data[nn]['xu'][n,0],traj.data[nn]['yu'][n,0],traj.data[nn]['zu'][n,0]))
