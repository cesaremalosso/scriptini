import numpy as np
import argparse
from subprocess import check_output, run
import glob as gl
import time
import multiprocessing as mp
BOHR = 0.529177249    # Bohr constant in Angstrom
TAU  = 0.5*4.8378e-5  # tau_CP constant in ps
HARTREE = 27.211386245988 #eV
eV=1.60217662 #10^-19 J
bar=1.0e-6/eV #eV/Ang^3

def   write_wanin(inputf,outf,pos,box,typedic,types ):
    nline = int(check_output(["wc", "-l", inputf]).decode("utf8").split()[0])
    with open(inputf,'r') as rf, open(outf,'w') as of :
        control=True
        for i in range(nline):
           line = rf.readline()
           lines = line.split()
           ls = len(lines)
           if(ls>1 and lines[0]=='begin' and lines[1]=='atoms_cart'):
             of.write(line)
             of.write('bohr \n')
             for i in range(len(pos)):
                s='{}  {:20.15f}  {:20.15f}  {:20.15f} \n'.format(typedic[types[i]],pos[i,0],pos[i,1],pos[i,2])
                of.write(s)
           elif(ls>1 and lines[0]=='begin' and lines[1]=='unit_cell_cart'):
             of.write(line)
             of.write('bohr \n')
             of.write('{:12.13f}  0.0  0.0 \n'.format(box[0]))
             of.write('0.0  {:12.13f}  0.0 \n'.format(box[0]))
             of.write('0.0  0.0  {:12.13f} \n'.format(box[0]))
           else:
             of.write(line)
    return 

def   write_pws(inputf,outf,pos,box,typedic,types ):
    nline = int(check_output(["wc", "-l", inputf]).decode("utf8").split()[0])
    with open(inputf,'r') as rf, open(outf,'w') as of :
        control=True
        for i in range(nline):
           line = rf.readline()
           lines = line.split()
           ls = len(lines)
           if(ls>0 and lines[0]=='ATOMIC_POSITIONS'):
             of.write(line)
             for i in range(len(pos)):
                s='{}  {:20.15f}  {:20.15f}  {:20.15f} \n'.format(typedic[types[i]],pos[i,0],pos[i,1],pos[i,2])
                of.write(s)
           elif(ls>0 and (lines[0]=='celldm(1)' or lines[0]=='celldm(1)=')):
             of.write('    celldm(1)= {:12.13f} \n'.format(box[0]))
           else:
             of.write(line)
    return 



def read_raw(folder,nstep,start):
    flist = [folder+'/coord.raw',folder+'/box.raw',folder+'type.raw']
    with open(folder+'/coord.raw','r') as rc,open(folder+'/box.raw') as rb,open(folder+'/type.raw') as rt : 
       tline = rt.readline().split()  
       natoms = len(tline)
       print('natoms in read_raw ', natoms)
       
       pos = np.zeros((nstep,natoms,3))
       box = np.zeros((nstep,9))
       #types 
       types = np.zeros(len(tline),dtype=int)
       for i in range(len(tline)):
          types[i] = int(tline[i])
       for i in range(len(tline)):
              types[i] = int(tline[i])

       #read box and positions
       #first skip start lines
       for i in range(start):
           cline = rc.readline()  
           bline = rb.readline()  
       #then read the following nstep lines
       for istep in range(nstep): 
           cline = rc.readline().split()  
           bline = rb.readline().split()  
           #pos
           for i in range(natoms):
              pos[istep,i,0] = float(cline[i*3+0])
              pos[istep,i,1] = float(cline[i*3+1])
              pos[istep,i,2] = float(cline[i*3+2])
           # box
           for i in range(len(bline)):
              box[istep,i] = float(bline[i])
           
    return pos , box,types

def gen_folder(outdir,istep,inputf,prefix,pos,box,typedic,types):
     outdirit=outdir+str(istep)
     if(len(gl.glob(outdirit)) == 0): run(["mkdir", "-p", outdirit+'/'])
     
     write_pws(inputf+'.scf',  outdirit+'/'+prefix+'.scf', pos[istep,:]*1./BOHR,box[istep,[0,4,8]]*1./BOHR,typedic , types )
#     write_pws(inputf+'.nscf', outdirit+'/'+prefix+'.nscf',pos[istep,:]*1./BOHR,box[istep,[0,4,8]]*1./BOHR,typedic , types )
#     write_wanin(inputf+'.win',outdirit+'/'+prefix+'.win', pos[istep,:]*1./BOHR,box[istep,[0,4,8]]*1./BOHR,typedic , types  )
     return

if( __name__ == "__main__"):
   parser = argparse.ArgumentParser(description = 'Convert .raw file in snapsho for pw')
   parser.add_argument('-d', '--directory',
   		type = str,
   		required = False,
   		help = 'Directory or the *raw files',
   		default = './')
   parser.add_argument('-p', '--prefix',
   		type = str,
   		required = True,
   		help = 'Prefix of the raw files.')
   parser.add_argument('--dout',
		type = str,
   		required = False,
   		help = 'Directory of the output file',
   		default = './')
   parser.add_argument('-i', '--input',
   		type = str,
   		required = True,
   		help = 'prefix for the input-examples of scf and nscf calculations. The files must be called prefix.scf prefix.nscf')
   parser.add_argument('-n', '--nstep',
   		type = int,
   		required = False,
                default = None,
   		help = 'number of steps')
   parser.add_argument( '--start',
   		type = int,
   		required = True,
                default = 0,
   		help = 'number of steps to skip at the begining of the file')
   parser.add_argument('-s', '--species',
		nargs = '*',
		type = str,
		required = True,
		help = 'Sequence of atomic species in the simulation (in the same order as in the types.raw input).')
   parser.add_argument( '--init_folder',
		type = int,
		required = False,
                default = 0,
		help = 'name of the initial folder.')
   parser.add_argument( '--parallel',
                action = 'store_true',
		required = False,
                default = False,
		help = 'enable parallel')
   parser.add_argument( '--npool',
		type = int,
		required = False,
                default = None,
		help = 'number of pools (processes) fo parallel execution, default total # of cpus.')
                
   args = parser.parse_args()
   
   directory = args.directory
   prefix = args.prefix
   outdir = args.dout
   inputf = args.input
   nstep = args.nstep
   init_folder = args.init_folder
   species = args.species
   parallel =args.parallel
   npool = args.npool
   start = args.start
  

   print("BEGIN")
   print("parallel= ",parallel)
   t1 = time.time()

   if(nstep is None):
     nstep = int(check_output(["wc", "-l", directory + '/box.raw']).decode("utf8").split()[0])
     nstep=nstep-start
   print('nstep',nstep) 
   if (npool is None):
     npool=mp.cpu_count()
   pos,box,types = read_raw(directory,nstep,start)
   typess =list(dict.fromkeys(types))
   typess.sort()
   typedic = {}
   for i,j in zip(typess,species):
     typedic[i] = j   
   t2=time.time()
   print("time to read {}".format(t2-t1))
   if(parallel):
      pool = mp.Pool(npool)
      for istep in range(nstep):  
         pool.apply_async(gen_folder,args=(outdir,istep+init_folder,inputf,prefix,pos,box,typedic,types))
         #gen_folder(outdir,istep+init_folder,inputf,prefix,pos,box,typedic,types)
      pool.close()
      pool.join()
   else:
      for istep in range(nstep):
         gen_folder(outdir,istep+init_folder,inputf,prefix,pos,box,typedic,types)

   t3=time.time()
   print('time to write {}'.format(t3-t2))
   print('total time {}'.format(t3-t1))
   print('END')


