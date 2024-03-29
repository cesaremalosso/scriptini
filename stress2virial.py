#!/bin/env python
##### coding=utf-8

import numpy as np
import argparse
import warnings
from subprocess import check_output
import re



mass = [1,16]
#mass[1] = 1
#mass[2] = 16
vol = 12.425508894894122**3
natom = 192
eV = 6.24150974#e18 #J2eV
NA = 6.02214076#e23
out_virial_raw = open('virial.raw','w')

typ = np.loadtxt('type.raw', dtype=int)
with open('stress.raw','r') as filestr, open('vel.raw', 'r') as filevel:
    for iline, line in enumerate(filestr):
        # linestr = filestr.readline()
        # convert from bar to eV
        values = np.array(line.split(), dtype=float) * vol# * 1e-7
        stress = np.reshape(values, (3, 3))

        linevel = filevel.readline()
        values = np.array(linevel.split(), dtype =float)
        vel = np.reshape(values, (natom, 3))

        kin = np.zeros((3,3))
        for iatom, tatom in enumerate(typ):
            kin += np.outer(vel[iatom], vel[iatom]) * mass[tatom] / NA * eV * 1e-4

        virial = stress + kin

        np.savetxt(out_virial_raw, np.reshape(virial, 9), newline=" ")
        out_virial_raw.write('\n')
out_virial_raw.close()
