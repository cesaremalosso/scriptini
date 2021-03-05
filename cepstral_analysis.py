#!/usr/bin/python

import sys
sys.path.append('/marconi/home/userexternal/cmalosso/thermocepstrum/themocepstrum')
import thermocepstrum as tc
import numpy as np
import argparse


# data in input = la temperatura
def block_analysis(data, tmax, dt_fs, vol):
    jen = data[:,3:]
    temp = data[:,2]

    Nstep = int(np.rint(tmax / (dt_fs * 1e-3)))
    maxrows = np.size(jen, 0)
    Ncurrs = maxrows / Nstep

    psds1 = 0
    psds2 = 0

    T = []

    for ij in range(int(Ncurrs)):
        init = Nstep * ij
        end = Nstep * (ij + 1) if Nstep * (ij + 1) < jen.shape[0] else jen.shape[0]

        Tmean = np.mean(temp[init:end])
        T.append(Tmean)

        jj = tc.HeatCurrent(j=jen[init:end], DT_FS=dt_fs, TEMPERATURE=Tmean, units='metal_vis', VOLUME=vol,
                            PSD_FILTER_W=0.05)

        # S(w=0) average
        psds1 += (jj.psd[0] - psds1) / (ij + 1)
        psds2 += (jj.psd[0] ** 2 - psds2) / (ij + 1)

    psdstd = np.sqrt((psds2 - psds1 ** 2) / (Ncurrs - 1))

    result = [T, psds1 * jj.kappa_scale / 2, psdstd * jj.kappa_scale / 2]
    return result

def main():
    parser = argparse.ArgumentParser(description='Block analysis script')
    parser.add_argument('--dt',
                        type=float,
                        required=True,
                        help='time step in femtoseconds')
    parser.add_argument('-v','--vol',
                        type=float,
                        required=True,
                        help='volume')
    parser.add_argument('--first',
                        type=int,
                        required=True,
                        help='first file')
    parser.add_argument('--last',
                        type=int,
                        required=True,
                        help='last file')
    parser.add_argument('--timelist',
                        type = int,
                        nargs='+',
                        required=True,
                        help='list of the size of blocks for block-analysis')
    args = parser.parse_args()

    first = args.first
    last = args.last

    lista = []
    for i in range(first,last+1):
        lista.append(np.loadtxt('stress{}.{}fs.out'.format(i,int(args.dt)),skiprows=1))
    data = np.concatenate(lista)


    for time in args.timelist:
        block_result = block_analysis(data, time, args.dt, args.vol)
        np.savetxt('{}ps_result.out'.format(time), block_result, fmt='%d')

if __name__ == "__main__":
    main()


