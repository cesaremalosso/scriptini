#!/bin/env python
##### coding=utf-8
import numpy as np
import argparse
import warnings
from subprocess import check_output
import re

BOHR = 0.529177249    # Bohr constant in Angstrom
#TAU  = 4.8378e-5  # tau_PW constant in ps
TAU  = 0.5*4.8378e-5  # tau_CP constant in ps
HARTREE = 27.211386245988 #eV
eV=1.60217662 #10^-19 J
bar=1.0e-6/eV #eV/Ang^3

def read_file_pos_vel(prefix, natoms, nstep=None):
	"""
	Legge i file di output di quantum espresso (cartella dove sono posizionati i vari restart)
	Per esempio, se il prefisso is KCl-512:
	namepos is KCl-512.pos
	namevel is KCl-512.vel
	nstep is il numero di timestep da leggere (NON is determinato automaticamente!)
	natoms is il numero di atomi nella simulazione (NON is determinato automaticamente!)
	
	ritorna una lista, una per ogni timestep, con il contenuto:
	    [timestep,tempo]       [ posizioni ]               [velocita']
	dove [timestep,tempo] is una coppia di numeri, posizioni is un array numpy con le posizioni,
	e velocita' is un array numpy
	"""	
	def file_length( filename ):
	  i = -1
	  blank = 0
	  with open(filename) as f:
	    for i, l in enumerate(f,1):
	      if len(l) == 1:
	        blank += 1
	      pass
	  return i - blank
	
	if nstep is None:
		nlines = int(check_output(["wc", "-l", prefix + '.pos']).decode("utf8").split()[0])
		nstep = nlines // (natoms+1)
		#nstep = file_length(prefix + '.evp') - 1
		print("nstep not set: using all the steps in .pos file: nstep = {}".format(nstep))

# ToDo: possibilita' di leggere i file dal fondo
#	if nstep < 0:
#		reverse = True
#		nstep = -nstep
#	else:
#		reverse = False
#
	data = {}
	
	filethe = open(prefix + '.evp', 'r')
	data['step']  = np.zeros(nstep, dtype = np.int64)
	data['time']  = np.zeros(nstep, dtype = np.float64)
	data['ekinc'] = np.zeros(nstep, dtype = np.float64)
	data['Tcell'] = np.zeros(nstep, dtype = np.float64)
	data['Tion']  = np.zeros(nstep, dtype = np.float64)
	data['econt'] = np.zeros(nstep, dtype = np.float64)
	data['epot'] = np.zeros(nstep, dtype = np.float64)
	#filethe.readline()  # skip first line
	
	filepos = open(prefix + '.pos', 'r')
	data['pos']  = np.zeros((nstep,natoms,3), dtype = np.float64)
	
	filevel = open(prefix + '.vel', 'r')
	data['vel']  = np.zeros((nstep,natoms,3), dtype = np.float64)
	
	try:
		filefor = open(prefix + '.for', 'r')
		read_force = True
		data['for']  = np.zeros((nstep,natoms,3), dtype = np.float64)
	except IOError as err:
		err = re.sub(r'\[.*\]', '', str(err))
		print('Warning!' + err + '. The .for file is not present: it will be ignored')
		read_force = False
	
	try:
		filestr = open(prefix + '.str', 'r')
		read_stress = True
		data['str']  = np.zeros((nstep,natoms,9), dtype = np.float64)
	except IOError as err:
		err = re.sub(r'\[.*\]', '', str(err))
		print('Warning!' + err + '. The .str file is not present: it will be ignored')
		read_stress = False
	
	filecel = open(prefix + '.cel', 'r')
	data['cell'] = np.zeros((nstep, 9), dtype = np.float64)
	
	istep = 0

	while (istep < nstep):
		linethe = filethe.readline()
		#print(linethe)
		linepos = filepos.readline()
		linevel = filevel.readline()
		if read_force: linefor = filefor.readline()
		if read_stress: linestr = filestr.readline()
		linecel = filecel.readline()
		if (len(linethe)==0) or (len(linepos)==0) or (len(linevel)==0) or (len(linecel)==0):  # EOF
			if read_force:
				if (len(linefor)==0):
					raise RuntimeError("End Of file")
			if read_stress:
				if (len(linestr)==0): raise RuntimeError("End Of file")
	
		# controllo per commenti 
		if (linethe.split()[0] == '#') :
			print("Comment found in {}.evp at line {}. Please check that this is correct.".format(prefix, istep+1))
			linethe = filethe.readline()
		# lettura thermo
		values = np.array(linethe.split(), dtype = np.float)
		if len(values):
			#print istep, values[0], len(data['step'])
			data['step'][istep]  = int(values[0])
			data['time'][istep]  = values[1]
			if istep == 1:
				deltat = data['time'][1]-data['time'][0] 
			data['ekinc'][istep] = values[2]
			data['Tcell'][istep] = values[3]
			data['Tion'][istep]  = values[4]
			data['epot'][istep]  = values[5]
			data['econt'][istep] = values[8]
		else:
			istep -= 1
	
		# lettura posizioni
		#values = np.array(linepos.split(), dtype = np.float)
		values = linepos.split()
		#print linepos
		#print values, data['step'][istep]
		if len(values):
			if (data['step'][istep] != int(values[0]) ):
				print(data['step'][istep], int(values[0]))
				raise RuntimeError("Different timesteps between files of positions and thermo")
			for iatom in range(natoms):
				linepos = filepos.readline()
				values = np.array(linepos.split())
				data['pos'][istep,iatom,:] = values[:]
	        
		#lettura velocity
		#values = np.array(linevel.split(), dtype=np.float)
		values = linevel.split()
		#print values,data[0][istep]
		if len(values):
			if (data['step'][istep] != int(values[0]) ):
				print(data['step'][istep], int(values[0]))
				raise RuntimeError("Different timesteps between files of velocity and thermo")
			for iatom in range(natoms):
				linevel = filevel.readline()
				values = np.array(linevel.split())
				data['vel'][istep,iatom,:] = values[:]
	    
		#lettura forza
		if read_force:
			values = linefor.split()
			#values = np.array(linefor.split(), dtype=np.float)
			#print values,data[0][istep]
			if len(values):
				if (data['step'][istep] != int(values[0]) ):
					print(data['step'][istep], int(values[0]))
					raise RuntimeError("Different timesteps between files of forces and thermo")
				for iatom in range(natoms):
					linefor = filefor.readline()
					values = np.array(linefor.split())
					data['for'][istep,iatom,:] = values[:]
	    
		#lettura stress
		if read_stress:
			#values = np.array(linestr.split(), dtype=np.float64)
			values = linestr.split()
			#print values,data[0][istep]
			if len(values):
				if (data['step'][istep] != int(values[0]) ):
					print(data['step'][istep], int(values[0]))
					raise RuntimeError("Different timesteps between files of stress and thermo")
				for iiline in range(3):
					linestr = filestr.readline()
					values = np.array(linestr.split())
					data['str'][istep,iatom,3*iiline:3*iiline+3] = values[:]
	    
		#lettura cella
		#values = np.array(linecel.split(), dtype=np.float64)
		values = linecel.split()
		#print values, data['step'][istep]
		if len(values):
			if (data['step'][istep] != int(values[0]) ):
				print(data['step'][istep], int(values[0]))
				raise RuntimeError("Different timesteps between files of cell and thermo")
			for i in range(3):
				values = np.array(filecel.readline().split())
				data['cell'][istep, 3*i] = values[0]
				data['cell'][istep, 3*i+1] = values[1]
				data['cell'][istep, 3*i+2] = values[2]
	
		istep += 1
	return data


def write_xyz(outfile, data, natoms_per_type, type_names=None, xyz = False, vel = False, charge = None, tskip = 1, vcm = False, raw = False, shuffle = False):
	"""
	Scrive un file nel formato lammpstrj (facilmente leggibile da vmd).
	cp.x nell'output separa gli atomi per tipi. Questa funzione assume lo stesso ordine.
	outfile is il nome del file da scrivere.
	data is il risultato della chiamata a read_file_pos_vel
	l is la dimensione della cella cubica scritta nell'output. """

	## Conversion factors
	conv_pos = BOHR
	conv_vel = BOHR/TAU
	conv_for = HARTREE/BOHR
	conv_energy = HARTREE
	conv_virial = bar

	## Put data in variables: improve readability
	POS = data['pos'] * conv_pos
	VEL = data['vel'] * conv_vel
	if 'for' in data: FOR = data['for'] * conv_for
	CELL = data['cell'] * conv_pos
	STEP = data['step']
	TEMP = data['Tion']
	EPOT = data['epot'] * conv_energy
	if 'str' in data: VIR = data['str'] * conv_virial

	out_file = open(outfile, "w")
	
	if xyz:
		out_xyz = open(prefix + '.xyz', 'w')
	
	if analisi:
		out_anal = open(prefix + '.analisi', 'w')
	
	if charge is not None:
		out_j = open(prefix + '.current', 'w')
		out_j.write('step   Temp   c_jion[1]   c_jion[2]   c_jion[3] # e*Ang/ps\n')
	
	if vcm:
		out_vcm = []
		for ityp, typ in enumerate(species):
			out_vcm.append(open(prefix + '.{}.vcm'.format(typ), 'w'))
			out_vcm[ityp].write('step   Temp   c_vcm{spec:s}[1]   c_vcm{spec:s}[2]   c_vcm{spec:s}[3] # Ang/ps\n'.format(spec=typ))

	if raw:
		out_box_raw = open('box.raw', 'w')
		out_coord_raw = open('coord.raw', 'ba')
		out_coord_raw.truncate(0)
		if 'for' in data: 
			out_force_raw = open('force.raw', 'ba')
			out_force_raw.truncate(0)
		out_energy_raw = open('energy.raw', 'w')
		if 'str' in data:
			out_virial_raw = open('virial.raw', 'ba')
			out_virial_raw.truncate(0)
	
	
	#out_file.write("This Text is going to out file\nLook at it and see\n")
	nsteps = POS.shape[0]
	natoms = POS.shape[1]
	if (natoms != sum(natoms_per_type)):
		raise ValueError('Sum of number of atoms per type does not match the total number of atoms.')
	if type_names is None:
		type_names = map(str, np.arange(1, len(natoms_per_type)+1))
	else:
		if (len(natoms_per_type) != len(type_names)):
			raise ValueError('Number of type_names not compatible with natoms_per_type.')
	

	if vel:
		np.savetxt(prefix + '.atmvel', np.reshape(VEL, (nsteps, 3*natoms)))
		
	# if needed, shuffle data for 'raw' files generation
	steps_raw = np.arange(0, nsteps, dtype = np.int)
	if shuffle:
		np.random.shuffle(steps_raw)

	for itimestep in range(0, nsteps, tskip):

		if xyz:
			out_xyz.write('{}\n\n'.format(natoms))
		
		if analisi:
			out_anal.write('{}\n'.format(natoms))
			out_anal.write('{} {}\n'.format(0, CELL[itimestep,0]))
			out_anal.write('{} {}\n'.format(0, CELL[itimestep,4]))
			out_anal.write('{} {}\n'.format(0, CELL[itimestep,8]))
	
		out_file.write("ITEM: TIMESTEP\n")
		out_file.write("{}\n".format(int(round(STEP[itimestep]))))
		out_file.write("ITEM: NUMBER OF ATOMS\n")
		out_file.write("{}\n".format(natoms))
		out_file.write('ITEM: BOX BOUNDS pp pp pp\n')
		out_file.write('{} {}\n'.format(0, CELL[itimestep,0]))
		out_file.write('{} {}\n'.format(0, CELL[itimestep,4]))
		out_file.write('{} {}\n'.format(0, CELL[itimestep,8]))
		out_file.write('ITEM: ATOMS id type x y z vx vy vz\n')
		cumnattype = np.cumsum(np.append(0,natoms_per_type))
		
		jion = np.zeros(3, dtype = np.float)

		if vcm:
			vcom = np.zeros((len(species), 3), dtype = np.float)
			vcom2 = np.zeros((len(species), 3), dtype = np.float)
		
		if raw:
			itimestep_raw = steps_raw[itimestep]
			# generate, once and for all, the type.raw file
			if itimestep == 0:
				with open('type.raw', 'w') as f:
					for attype, nattype in enumerate(natoms_per_type):
						firstat = cumnattype[attype]
						lastat  = cumnattype[attype+1]
						for i, idat in enumerate(range(firstat, lastat)):
							f.write('{} '.format(attype))
			#out_coord_raw.write(np.reshape(POS[itimestep_raw, :, :], 3*natoms)) 
			#np.savetxt(out_coord_raw, np.reshape(POS[itimestep_raw, :, :], 3*natoms), newline = " ") 
			tosave = POS[itimestep_raw, :, :]
			sides = np.diag(np.reshape(CELL[itimestep_raw], (3, 3)))
			tosave = np.reshape(tosave%sides, 3*natoms) # TODO: implement PBC in more general cases
			np.savetxt(out_coord_raw, tosave, newline = " ")
			out_coord_raw.write('\n'.encode("utf-8"))
			#if 'for' in data: out_force_raw.write(np.reshape(FOR[itimestep_raw, :, :], 3*natoms)) 
			if 'for' in data:
				np.savetxt(out_force_raw, np.reshape(FOR[itimestep_raw, :, :], 3*natoms), newline = " ")
				out_force_raw.write('\n'.encode("utf-8"))
			if 'str' in data:
				np.savetxt(out_virial_raw, np.reshape(VIR[itimestep_raw, :, :], 9*natoms), newline = " ")
				out_virial_raw.write('\n'.encode("utf-8"))
			#out_box_raw.write(np.reshape(CELL[itimestep_raw, :, :], 9*natoms))
			for ibox in range(8):
				out_box_raw.write('{} '.format(CELL[itimestep_raw, ibox]))
			out_box_raw.write('{}\n'.format(CELL[itimestep_raw, -1]))
			out_energy_raw.write('{}\n'.format(EPOT[itimestep_raw]))

		for attype, nattype in enumerate(natoms_per_type):
			firstat = cumnattype[attype]
			lastat  = cumnattype[attype+1]
			for i, idat in enumerate(range(firstat, lastat)):
				out_file.write('{} {} {} {} {} {} {} {}\n'.format(idat+1, type_names[attype], \
					POS[itimestep,idat,0],
					POS[itimestep,idat,1],
					POS[itimestep,idat,2],
					VEL[itimestep,idat,0],
					VEL[itimestep,idat,1],
					VEL[itimestep,idat,2] \
					))
				if xyz:
					if vel:
						out_xyz.write('{:s} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f}\n'.format(
							type_names[attype],
							POS[itimestep,idat,0], 
							POS[itimestep,idat,1],
							POS[itimestep,idat,2],
							VEL[itimestep,idat,0],
							VEL[itimestep,idat,1], 
							VEL[itimestep,idat,2] \
							))
					else:
						out_xyz.write('{:s} {:15.10f} {:15.10f} {:15.10f}\n'.format(
							type_names[attype],
							POS[itimestep,idat,0], 
							POS[itimestep,idat,1],
							POS[itimestep,idat,2] \
							))
				if analisi:
					out_anal.write('{:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f}\n'.format(idat+1, attype, \
						POS[itimestep,idat,0],
						POS[itimestep,idat,1],
						POS[itimestep,idat,2],
						VEL[itimestep,idat,0],
						VEL[itimestep,idat,1],
						VEL[itimestep,idat,2] \
						))

				if charge is not None:

					jion += charge[attype] * VEL[itimestep, idat, :]

				if vcm:
					vcom[attype, :]  += VEL[itimestep, idat, :]
					vcom2[attype, :] += VEL[itimestep, idat, :]**2

		if vcm:
			for typ, typname in enumerate(species):
				out_vcm[typ].write('{:d} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f}\n'.format(int(round(STEP[itimestep])),
					TEMP[itimestep],
					vcom[typ][0],
					vcom[typ][1],
					vcom[typ][2],
					vcom2[typ][0],
					vcom2[typ][1],
					vcom2[typ][2] \
					))

		if charge is not None:
			out_j.write('{:d} {:15.10f} {:15.10f} {:15.10f} {:15.10f}\n'.format(int(round(STEP[itimestep])),
				TEMP[itimestep],
				jion[0],
				jion[1],
				jion[2] \
				)) 
	out_file.close()
	return
 
class lammpsfile:
    def __init__(self,natoms,prefix='dump',masses=None,ntypes=None,types=None,nbonds=None,nangles=None):
        self.prefix = prefix
        self.natoms = natoms
        self.nbonds = nbonds
        self.nangles = nangles
        self.ntypes = ntypes
        self.box = None
        self.pos= None
        self.masses = masses
        self.types = types
        if(self.natoms!=len(self.types)):
            print("self.natoms!=len(self.types) {} != {}".format(self.natoms,len(self.types)))
        self.step=0
        return
    def open_write(self):
        self.ftraj = open('./'+self.prefix+'.lammpstrj','w')
        self.fdata = open('./'+self.prefix+'.data','w') 
        return 
    def open_read(self):
        self.ftraj = open('./'+self.prefix+'.lammpstrj','r')
        self.fdata = open('./'+self.prefix+'.data','r') 
        return 
    def write_data_header(self):
        s= ' LAMMPS Description  \n \n'
        s+='     {} atoms\n'.format(self.natoms)
        s+='     {} bonds\n'.format(self.nbonds)
        s+='     {} angles\n \n'.format(self.nangles)
        s+='           {} atom types\n'.format(self.ntypes)
        s+='           0 bond types\n'
        s+='           0 angle types\n \n'
        s+='  {:2.1f}       {:10.9f}     xlo xhi \n'.format(0.0,self.box[0,0])
        s+='  {:2.1f}       {:10.9f}     ylo yhi \n'.format(0.0,self.box[1,1])
        s+='  {:2.1f}       {:10.9f}     zlo zhi \n \n'.format(0.0,self.box[2,2])
        s+=' Masses \n \n'
        #print(ntypes)
        for i in range(self.ntypes):
            s+='   {:5d}  {:10.9f} \n'.format(i+1,self.masses[i])
        s+=' Atoms\n \n' #full id molid type charge x y z
        self.fdata.write(s)
        return
        
        
    def close(self):
        self.fdata.close()
        self.ftraj.close()
        return
        
    def write_trj_header(self):
        s="ITEM: TIMESTEP\n"
        s+="{}\n".format(self.step)
        s+="ITEM: NUMBER OF ATOMS\n"
        s+="{}\n".format(self.natoms)
        s+="ITEM: BOX BOUNDS pp pp pp\n"
        s+="{} {}\n".format(0, self.box[0,0])
        s+="{} {}\n".format(0, self.box[1,1])
        s+="{} {}\n".format(0, self.box[2,2])
        s+="ITEM: ATOMS id type x y z\n"
        self.ftraj.write(s)
        return
    def write_step(self):
        for i in range(self.natoms):
            self.ftraj.write('{} {} {} {} {} \n'.format(i+1, int(self.types[i])+1, \
        self.pos[i,0],self.pos[i,1],self.pos[i,2]))
        return

    def read_step(self):
        line = self.ftraj.readline().split()
        line = self.ftraj.readline().split()
        self.step = float(line[0])
        line = self.ftraj.readline().split()
        line = self.ftraj.readline().split()
        self.natoms = int(line[0])
        line = self.ftraj.readline().split()
        self.box = np.zeros((3,3))
        self.pos = np.zeros((self.natoms,3))
        self.force = np.zeros((self.natoms,3))
        self.types = np.zeros(self.natoms,dtype=int)
        for i in range(3):
            line = self.ftraj.readline().split()
            self.box[i,i] = float(line[1])-float(line[0])
        line = self.ftraj.readline().split() 
        items = {}
        for i in range(2,len(line)):
            if(line[i]=='id'):
               items[line[i]] = i-2
            elif(line[i]=='type'):
               items[line[i]] = i-2
            elif(line[i][0]=='x'):
               items['x'] = i-2
            elif(line[i][0]=='y'):
               items['y'] = i-2
            elif(line[i][0]=='z'):
               items['z'] = i-2
            elif(line[i]=='fx'):
               items[line[i]] = i-2
            elif(line[i]=='fy'):
               items[line[i]] = i-2
            elif(line[i]=='fz'):
               items[line[i]] = i-2
               
        for i in range(self.natoms):
            line = self.ftraj.readline().split() 
            idx = int(line[items['id']])
            #print(i,idx)
            self.types[idx-1]   = int(line[items['type']])
            self.force[idx-1,0] = float(line[items['fx']])  
            self.force[idx-1,1] = float(line[items['fy']])  
            self.force[idx-1,2] = float(line[items['fz']])  
            self.pos[idx-1,0]   = float(line[items['x'] ]) 
            self.pos[idx-1,1]   = float(line[items['y'] ]) 
            self.pos[idx-1,2]   = float(line[items['z'] ]) 
            
        
        return
        
    def convert_box(self,bbox):
        self.box = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                self.box[i,j]=bbox[i*3+j]
        return
        
        
        
        
    
 
class rawfile:
    def __init__(self,ifforce,folder):
        self.rc = open(folder+'/coord.raw','r')
        self.rb = open(folder+'/box.raw'  ,'r')
        self.box = np.zeros(9)
        
        self.types = np.loadtxt(folder+'/type.raw')
        self.natoms = len(self.types)
        mini = np.amin(self.types)
        maxa = np.amax(self.types)
        self.ntype = int(maxa - mini +1)
        self.natoms_per_types = np.zeros(self.ntype)
        for i in self.types:
            self.natoms_per_types[int(i-mini)] += 1
        
        self.pos = np.zeros((self.natoms,3))
        self.rf=None
        if(ifforce):
           self.rf=open(folder+'/force.raw' ,'r')
           self.forc = np.zeros((self.natoms,3))
        return
    
    def read_step(self):
        cline = self.rc.readline().split()
        if(len(cline)==0):
            print("len(coordline)==0")
            return
        bline = self.rb.readline().split()
        if(len(bline)==0):
            print("len(boxline)==0")
            return
        if(self.rf is not None):
            fline = self.rf.readline().split()
            if(len(fline)==0):
                print("len(forceline)==0")
                return
        # pos, force
        for i in range(self.natoms):
            self.pos[i,0] = float(cline[i*3+0])
            self.pos[i,1] = float(cline[i*3+1])
            self.pos[i,2] = float(cline[i*3+2])
            if(self.rf is not None):
                self.forc[i,0] = float(fline[i*3+0])
                self.forc[i,1] = float(fline[i*3+1])
                self.forc[i,2] = float(fline[i*3+2])
        #box
        assert (len(bline) == 9),'len(bline) != 9 : {} != 9'.format(len(bline))
        for i in range(9):
            self.box[i] = float(bline[i])
        return
        
    def close(self):
        self.rc.close()
        self.rb.close()
        if(self.rf is not None): self.rf.close()
        return
    def __repr__(self):
       s='natoms = {}\n'.format(self.natoms)
       s+='ntype = {}\n'.format(self.ntype)

###########################################################################################################################################################################
### Parser
if( __name__ == "__main__"):
    parser = argparse.ArgumentParser(description = 'Convert raw     files into LAMMPS traj or data. The units are metal.')
    parser.add_argument('-d', '--directory',
    		type = str,
    		required = False,
    		help = 'Directory with the .raw files.',
    		default = './tmp')
    parser.add_argument('-p', '--prefix',
    		type = str,
    		required = False,
    		help = 'Prefix of the lammps filename.',
            default = 'lammps')
    parser.add_argument('-c', '--charge',
    		nargs = '*',
    		type = float,
    		required = False,
    		help = 'Oxidation number per species (in the same   order as in the ATOMIC_SPECIES card in the cp.x     input).')
    parser.add_argument('--nstep',
    		type = int,
    		default = 100,
    		help = 'Number of steps to convert.')
    parser.add_argument('--tskip',
    		type = int,
    		default = 1,
    		help = 'Write 1 every tskip steps.')
    parser.add_argument('--data',
    		help = 'Write the prefix.data file.',
    		action = 'store_true',
    		required = False,
    		default = False)
    parser.add_argument('--lmptraj',
            help = 'Write the prefix.lammpstrj file.',
            action = 'store_true',
            required = False,
            default = True)
    parser.add_argument('--force',
    		help = 'Write also the velocities in the .lammpstrj     file.',
    		action = 'store_true',
    		required = False,
    		default = False)
    parser.add_argument('--analisi',
    		help = 'Output data in analisi format.',
    		action = 'store_true',
    		default = False)
    
    
    args = parser.parse_args()
    
    folder = args.directory
    prefix = args.prefix
    #species = args.species
    #natm = args.natm
    nstep = args.nstep
    force = args.force
    analisi = args.analisi
    charge = args.charge
    tskip = args.tskip
    
    
    flist=[folder+'/coord.raw',folder+'/box.raw',folder+'type.raw']
    if(force):
        flist=[folder+'/coord.raw',folder+'/box.raw',folder+'type.raw',folder+'/force.raw']
    
    print('Nstep = {}'.format(nstep))
        
    
    print('Reading {}/*raw files ...'.format(folder, prefix))
    raw = rawfile(ifforce=force,folder=folder)
    lammps = lammpsfile(prefix=prefix,natoms=raw.natoms,ntypes=raw.ntype,types=raw.types)
    lammps.open_write() 
    for i in range(nstep):
        raw.read_step()
        lammps.pos=np.copy(raw.pos)
        lammps.convert_box(bbox=raw.box)
        lammps.write_trj_header()
        lammps.write_step()
        lammps.step += 1
    lammps.close()
    raw.close()
    print('Done.')
