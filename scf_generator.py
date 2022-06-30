import numpy as np
import argparse

BOHR = 0.529177249    # Bohr constant in Angstrom



def write_input():
    with open(calculation+'.in', 'w') as fileout:
        fileout.writelines('&CONTROL, \n')
        fileout.writelines('title = "{}", \n'.format(prefix))
        fileout.writelines('restart_mode = "from_scratch", \n')
        fileout.writelines('verbosity = "medium"", \n')
        fileout.writelines('prefix = "{}", \n'.format(prefix))
        fileout.writelines('wf_collect = .false., \n')
        fileout.writelines('out_dir = "./out", \n')
        fileout.writelines('pseudo_dur = "{}", \n'.format(pseudo_dir))
        fileout.writelines('tprnfor = .true., \n')
        fileout.writelines('tstress = .true., \n')
        fileout.writelines('max_seconds = 86000, \n')
        fileout.writelines('etot_conv_thr = {}", \n'.format(econv))
        if calculation == 'vc-relax':
            fileout.writelines('forc_conv_thr = {}", \n'.format(fconv))
        fileout.writelines('/ \n')
        fileout.writelines('&SYSTEM \n')
        fileout.writelines('ibrav = {}, \n'.format(ibrav))
        fileout.writelines('nat = {}, \n'.format(natoms))
        fileout.writelines('ntyp = {}, \n'.format(ntypes))
        fileout.writelines('ecutwfc = {}, \n'.format(wfccut))
        fileout.writelines('ecutrho = {}, \n'.format(rhocut))
        fileout.writelines('ntyp = {}, \n'.format(ntypes))
        fileout.writelines('ntyp = {}, \n'.format(ntypes))
        fileout.writelines('ntyp = {}, \n'.format(ntypes))


###########################################################################################################################################################################
### Parser

parser = argparse.ArgumentParser(description = 'Generate a scf.in input file')
parser.add_argument('-d', '--directory',
                    type = str,
                    required = False,
                    help = 'Target directory',
                    default = './')
parser.add_argument('-p', '--prefix',
                    type = str,
                    required = False,
                    help = 'Prefix of the filename.',
                    default = 'scf')
parser.add_argument('-s', '--species',
                    nargs = '*',
                    type = str,
                    required = True,
                    help = 'Sequence of atomic species in the simulation.')
parser.add_argument('-n', '--natm',
                    type = int,
                    required = True,
                    help = 'Number of atoms per species.')
parser.add_argument('-c', '--charge',
                    nargs = '*',
                    type = float,
                    required = True,
                    help = 'Oxidation number per species (in the same order as in -s --species card).')
parser.add_argument('--calculation',
                    help = 'Kind of calculation you want to perform.',
                    required = True)
parser.add_argument('--pseudo',
                    nargs = '*',
                    type = str,
                    required = True,
                    help = 'Name of the pseudo files (in the same order as in -s --species card)')
parser.add_argument('--pseudo_dir',
                    type = str,
                    required = True,
                    help = 'Directory containing the pseudo files')

