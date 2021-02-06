#!/bin/env python

import numpy as np
import sys

name_input=sys.argv[1]
name_output=sys.argv[2]
number_particles=sys.argv[3]

out=open(str(name_output),'w')

out.write('Lammps input Description \n')
out.write('\n')

fp=open(str(name_input))
for i,line in enumerate(fp):
    if i==2:
        out.write(line)
    elif i==3:
        out.write(line)
    elif i==9:
        out.write('\n')
        out.write(line)
    elif 9<i<17:
        out.write(line)
    elif i==32:
        out.write('\n')
        out.write('Atoms # atomic \n')
        out.write('\n')
    elif 32<i< int(number_particles)+33:
        a=line.split()
        out.write(str(a[0]+'  '+str(a[2])+'  '+str(a[4])+'  '+str(a[5])+'  '+str(a[6]))+'\n')
out.close()
