# coding=utf-8

import numpy as np
import pandas as pd
import sys
def wrap_coord(lb,x):
    while(x > lb or x < 0):
        x = x - np.sign(x) *lb
    return

def read_write_dat(infile,outdir='./',step=-1):
    """
      INPUT: 
            infile = input file name
            outfile = output file name
            step = max number of timestep (if negative use all timesteps)
      OUTPUT:
    """
    snap = {}

    with open(infile, "r") as in_file, open(outdir + 'box.raw','w') as boxf, open(outdir + 'coord.raw','w') as coordf ,\
        open(outdir + 'force.raw','w') as forcef, open(outdir + 'Nstep.data','w') as stepf :
         end = True
         line = in_file.readline()
         if(line == ''):
             print('first line empty')
             return
         dt=0
         while (end):
                 dt += 1
                 timestep = int(in_file.readline().split()[0])
                 stepf.write('{} \n'.format(timestep))
                 line = in_file.readline()
                 natoms=int(in_file.readline().split()[0])
                 snap['coord'] = np.zeros((natoms,3))
                 snap['force'] = np.zeros((natoms,3))
                 snap['type'] = np.zeros(natoms,dtype=int)
                 old = np.zeros(natoms)
                 line = in_file.readline()
                 box = np.zeros(3)
                 llbox = np.zeros(3)
                 for i in range(3):
                    line = in_file.readline().split()
                    llbox[i] = float(line[0])
                    box[i] = float(line[1]) - llbox[i]
                 boxf.write('{} 0.0 0.0 0.0 {} 0.0 0.0 0.0 {} \n'.format(box[0],box[1],box[2]))
                 llbox = -llbox
                 line = in_file.readline().split()[2:]
                 #TODO understand what is dumped
                 for iatom in np.arange(natoms):
                      line =  in_file.readline().split()
                      snap['coord'][iatom,0] =  float(line[2]) + llbox[0]
                      snap['coord'][iatom,1] =  float(line[3]) + llbox[1]
                      snap['coord'][iatom,2] =  float(line[4]) + llbox[2]
                      for ii in range(3):
                          wrap_coord(box[ii],snap['coord'][iatom,ii])
                      snap['force'][iatom,0] =  float(line[8])
                      snap['force'][iatom,1] =  float(line[9])
                      snap['force'][iatom,2] =  float(line[10])
                      snap['type'][iatom]  =  int(line[1])-1
                      coordf.write('{} {} {} '.format(snap['coord'][iatom,0],snap['coord'][iatom,1],snap['coord'][iatom,2]))
                      forcef.write('{} {} {} '.format(snap['force'][iatom,0],snap['force'][iatom,1],snap['force'][iatom,2]))
                     # if (old[iatom] != snap['type'][iatom] and dt > 1 ):
                     #     print('old[',iatom,'] != snap[\'type\'] [',iatom,']  ->  ',old[iatom], ' != ',  snap['type'] [iatom])
                     #     print('at step ',dt)
                     #     return
                 old[:] = snap['type'][:]
                 coordf.write('\n')
                 forcef.write('\n')
                 line = in_file.readline()
                 if(line == '' or (step == dt and step >0)):
                     end=False

    return


def write_ener(enerfile,outdir='./', stepf='./Nstep.data'):
    ener = np.loadtxt(enerfile,skiprows=1)  
    np.savetxt(outdir + 'energy.raw', ener[:,2])
    #ener = pd.read_csv(enerfile, sep=None, engine='python').to_dict()
    #en = dict(zip(int(ener['Step']),ener['PotEng']))

    #steps = pd.read_csv(stepf,sep=None, engine='python').to_numpy(dtype='int')
    #with open(outdir + 'energy.raw') as efile:
    #    for i in range(steps.shape[0]):
    #       efile.write('{} \n'.format(en[i]))


    return

if len(sys.argv)<3:
    print ("Uso: {} [prefix_file_in] [file_out] [enerfile] [stop (optional)] ".format(sys.argv[0]))
    print ("script to read a lammps dump file and write them on a .raw file without any other line")
    print ("the structure of the dump is ITEM: ATOMS id type xu yu zu vx vy vz fx fy fz")
    exit(-1)

stop=-1
if (len(sys.argv)==5):
   stop=int(sys.argv[4])

read_write_dat(sys.argv[1],sys.argv[2],stop)
#write_ener(enerfile=sys.argv[3],outdir=sys.argv[2],stepf=sys.argv[2] + '/Nstep.data')

