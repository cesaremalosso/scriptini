import numpy as np
import argparse

def read_pw(infile):
    with open(infile, 'r') as inf:
        dd = True
        step = 0
        while (dd):
            dd+=1
            ll = inf.readline().split()
            leng = len(ll)
            if (leng > 3 and ll[0] == 'total' and ll[1] == 'stress'):
                stress = np.zeros((3, 3))
                for i in range(3):
                    ll = inf.readline().split()
                    stress[i, 0] = ll[3]
                    stress[i, 1] = ll[4]
                    stress[i, 2] = ll[5]
                dd = False
            else:
                if (step >= 3000):
                    fff = open('errore', 'w+')
                    fff.write('error  in pw')
                    fff.close()
                    dd = False

    return stress

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

    stress = read_pw(infile)

    vfile = outdir + '/stress_pw.out'
    with open(vfile, 'w') as vf:
        for i in range(3):
            vf.write(' {:18.15e}  {:18.15e}  {:18.15e} \n'.format(stress[i, 0], stress[i, 1], stress[i, 2]))
        vf.write('\n')
    return


if (__name__ == "__main__"):
    main()