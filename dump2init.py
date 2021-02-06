import sys
sys.path.append("/home/cmalosso/scriptini")
import read_lammps_dump as rd

a_B = 0.529177210903
name_input=sys.argv[1]
name_output=sys.argv[2]
#units = sys.argv[3]
every_n=int(sys.argv[3])
h_type=int(sys.argv[4])
o_type=int(sys.argv[5])

traj = rd.LAMMPS_Dump(name_input)

number_of_POSCARS=traj.TOT_TIMESTEPS//every_n
traj.read_timesteps((traj.FIRST_TIMESTEP,traj.LAST_TIMESTEP+traj.DELTA_TIMESTEP,traj.DELTA_TIMESTEP*every_n))

for nn in range(number_of_POSCARS):
	out=open(str(name_output)+str(nn),'w')
	
	out.write('POSCAR generated with dump2init \n')
	out.write('1.0 \n')
#	if (units == 'bohr'):
#		out.write(str(a_B)+'\n')
#	else:
#		print('only bohr to metal unit implemented')
	out.write('{}   {:.9f}   {:.9f}\n'.format(traj.BOX_BOUNDS[0,1],0,0))
	out.write('{:.9f}   {}   {:.9f}\n'.format(0,traj.BOX_BOUNDS[0,1],0))
	out.write('{:.9f}   {:.9f}   {}\n'.format(0,0,traj.BOX_BOUNDS[0,1]))
	
	out.write('H      O\n')
	
	out.write('{}   {}\n'.format(traj.NATOMS*2//3,traj.NATOMS//3))
	
	out.write('Cartesian\n')
	
	for n in range(traj.NATOMS):
		if (traj.data[nn]['type'][n]==h_type):
			out.write('{}   {}   {}\n'.format(traj.data[nn]['x'][n,0],traj.data[nn]['y'][n,0],traj.data[nn]['z'][n,0]))
	for n in range(traj.NATOMS):
		if (traj.data[nn]['type'][n]==o_type):
			out.write('{}   {}   {}\n'.format(traj.data[nn]['x'][n,0],traj.data[nn]['y'][n,0],traj.data[nn]['z'][n,0]))
	
	out.close()
