

import numpy as np 
import sys
import subprocess as sb
import os

Ha=27.2114 #eV
Ry= 13.6056980659 #eV
rBohr=0.529177249 #Ang 
eV=1.60217662 #10^-19 J
bar=1.0e-6/eV #eV/Ang^3
force_c = Ry/rBohr

def main(inputf,outputf,coeff_list,natoms=375,debug=False):
   data = np.loadtxt(inputf,comments='#')
   if(debug): print(data)
   lr = data[:,-1]
   w_f = np.zeros(lr.shape)
   w_e = np.zeros(lr.shape)
   w_e = np.zeros(lr.shape)
   start_lr = lr[0]
   batch = data[:,0]
   start_pref_e=float(coeff_list[0])
   limit_pref_e=float(coeff_list[1])
   start_pref_f=float(coeff_list[2])
   limit_pref_f=float(coeff_list[3])
   start_pref_v=float(coeff_list[4])
   limit_pref_v=float(coeff_list[5])

   w_f = start_pref_f * ( lr / start_lr ) + limit_pref_f * ( 1 - lr / start_lr ) 
   w_e = start_pref_e * ( lr / start_lr ) + limit_pref_e * ( 1 - lr / start_lr ) 
   w_v = start_pref_v * ( lr / start_lr ) + limit_pref_v * ( 1 - lr / start_lr ) 
   data[:,3] = data[:,3]*data[:,3]*natoms*w_e
   data[:,4] = data[:,4]*data[:,4]*natoms*w_e
   data[:,5] = data[:,5]*data[:,5]*w_f
   data[:,6] = data[:,6]*data[:,6]*w_f
   np.savetxt(outputf,np.column_stack((data,w_e,w_f,w_v)),fmt=('%3.5e'),header='batch      l2_tst    l2_trn    l2_e_tst  l2_e_trn    l2_f_tst  l2_f_trn         lr')
   return w_f,w_e,w_v
   


if (__name__ == '__main__'):
  import argparse
  parser = argparse.ArgumentParser(description=' compute the coefficient of the training')
  parser.add_argument('-i',dest='input', metavar='infile', type=str, \
                       help='file with step to read',default='lcurve.out')
  parser.add_argument('-o',dest='output', metavar='infile', type=str,\
                    help='file to write failed steps',default='coeff_lcurve.out')
  parser.add_argument('-c','--coeff', dest='lists',nargs='+', help='coefficients list [pe_in pe_fin pf_in pf_fin pv_in pv_fin]', default=[0.,0.,0.,0.,0.,0.])
  parser.add_argument('-d', dest='debug', metavar='db',type=bool,\
                    help='debug mode', default=False)
  parser.add_argument('-n', dest='natoms', metavar='na',type=int,\
                    help='number of atoms (Default: 375)', default=375)
  args = parser.parse_args()
  #print(args.input)
  coeff_list=args.lists  
  if(len(coeff_list)<6):
     for i in range(6-len(coeff_list)):
        coeff_list.append(0.)

  w_f,w_e,w_v=main(inputf=args.input,outputf=args.output, coeff_list=coeff_list,natoms=args.natoms,debug=args.debug)
   
