import numpy as np
import argparse

def wrap_coord(lpart,lbox):
    periodic=True
    while(periodic) :
      if(lpart < 0):
          lpart = lpart + lbox
      elif (lpart > lbox):
          lpart = lpart - lbox
      else:
          periodic=False
    return lpart


def read_ewald(inforce='energy_force_ewald.dat'):
  try: 
   with open(inforce,'r') as inf:
      ll=inf.readline().split()
      energy = float(ll[0])
   force = np.loadtxt(inforce,skiprows=1)
  except FileNotFoundError:
   fff = open('errore','w+')
   fff.write('errore ewald no file')
   fff.close()
   energy = None
   force = np.zeros((4,4))
   
   
  return energy, force[:,1:]

def read_pw(infile,logfile='log.log'):
   atom=None
   box=None
   with open(infile,'r') as inf:
     dd=True
     step=0
     while(dd):
       step+=1
       ll=inf.readline().split()
       leng=len(ll)
       if(leng>2 and ll[2]=='atoms/cell'):   # number of atoms/cell
          nat=int(ll[4])
       elif(leng>3 and ll[0]=='!'):   # 
          print('energy')
          energy= float(ll[4])
       elif (leng>3 and ll[0]=='celldm(1)=' ):
          print('box')
          box = np.diag(np.ones(3)*float(ll[1])).reshape(9)
       elif (leng>3 and ll[0]=='site' and ll[3]=='positions'):
          print(ll)
          atom = np.zeros((nat,3))
          for i in range(nat):
             ll=inf.readline().split()
             iatom = int(ll[0])-1
             atom[iatom,0] = float(ll[6])
             atom[iatom,1] = float(ll[7])
             atom[iatom,2] = float(ll[8])
          
       elif(leng>3 and ll[0]=='Forces' and ll[1]=='acting'):
          print(ll)
          ll=inf.readline().split()
          force=np.zeros((nat,4))
          for i in range(nat):
             ll=inf.readline().split()
             iatom = int(ll[1])-1
             force[iatom,0] = float(ll[6])
             force[iatom,1] = float(ll[7])
             force[iatom,2] = float(ll[8])
             force[iatom,3] = float(ll[3])*1.
          dd=False
       else: 
         if(step>=1000):
           fff = open('errore','w+')
           fff.write('error  in pw')
           fff.close()
           dd=False
           nat = None
           energy = None
           force = None

     return nat , energy , force , atom , box
def  main():
   parser = argparse.ArgumentParser(description = 'write force.raw and energy.raw with ewald contributions')
   parser.add_argument('-i','--pwout',
                  type = str,
                  required = False,
                  help = 'output pw/qe file',
                   default = './scf.out')
   parser.add_argument( '--outdir',
               type = str,
               required = False,
               default = './',
               help = 'output folder where to put the .raw files')
   
   args = parser.parse_args()
   infile= args.pwout
   outdir = args.outdir

   nat , energy , force , atom , box = read_pw(infile)
   if (nat == None or energy == None ):
           fff = open('errore','w+')
           fff.write('error  in pw')
           fff.close()
           return
       
      
   
   
   efile = outdir + '/energy_pw.raw' 
   ffile = outdir + '/force_pw.raw'
   cfile = outdir + '/coord_pw.raw'
   bfile = outdir + '/box_pw.raw'
   
   convEne = 13.605693122994 # Ry to eV
   bohr_radius = 0.529177210903 # bohr to Ams
   convForce = convEne / bohr_radius 
   energy *= convEne
   force[:,:3] *= convForce
   box *= bohr_radius
   atom *= box[0]

   np.savetxt(bfile,box.reshape(1,-1))
   np.savetxt('type_pw.raw',force[:,3].reshape(1,-1)-1,fmt='%d')

   with open(efile,'w') as ef , open(ffile,'w') as ff, open(cfile,'w') as cf :
     ef.write('{} \n'.format(energy))
     for iat in range(nat) : 
        ff.write(' {:18.15e}  {:18.15e}  {:18.15e} '.format(force[iat,0],force[iat,1],force[iat,2]))
        cf.write(' {:18.15e}  {:18.15e}  {:18.15e} '.format(atom[iat,0],atom[iat,1],atom[iat,2]))
     ff.write('\n')
     cf.write('\n')
   
   return 

if( __name__ == "__main__"):
  main()
