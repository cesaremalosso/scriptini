import numpy as np
import argparse

BOHR = 0.529177249e-10 # Bohr constant in Metres
Ry = 2.1798723611035e-18 # Rydberg constan in Joules

tobar = Ry/BOHR**3/1e5

def read_pw(infile):
    with open(infile, 'r') as inf:
        dd = True
        step = 0
        while (dd):
            dd+=1
            ll = inf.readline().split()
            leng = len(ll)
            if (leng > 3 and ll[0] == 'total' and ll[1] == 'stress'):
#                pressure = ll[5]
                stress = np.zeros((3, 3))
                for i in range(3):
                    ll = inf.readline().split()
                    stress[i, 0] = float(ll[0])*tobar
                    stress[i, 1] = float(ll[1])*tobar
                    stress[i, 2] = float(ll[2])*tobar
                pressure = np.trace(stress)/3
                dd = False
            else:
                if (step >= 3000):
                    fff = open('errore', 'w+')
                    fff.write('error  in pw')
                    fff.close()
                    dd = False

    return stress, pressure

def main():
    parser = argparse.ArgumentParser(description='grep the stress in kbar from pws files')
    parser.add_argument('-i', '--pwout',
                        type=str,
                        required=False,
                        help='output pw/qe file',
                        default='./scf.out')
    parser.add_argument('--outdir',
                        type=str,
                        required=False,
                        default='./',
                        help='output folder where to put the stress file')
    args = parser.parse_args()
    infile = args.pwout
    outdir = args.outdir

    stress, pressure = read_pw(infile)

    vfile = outdir + '/stress_pw.out'
    pfile = outdir + '/pressure_pw.out'

    with open(pfile, 'w') as pf:
        pf.write('{} \n'.format(pressure))

    with open(vfile, 'w') as vf:
        for i in range(3):
            vf.write(' {:18.15e}  {:18.15e}  {:18.15e} \n'.format(stress[i, 0], stress[i, 1], stress[i, 2]))
        vf.write('\n')
    return


if (__name__ == "__main__"):
    main()
