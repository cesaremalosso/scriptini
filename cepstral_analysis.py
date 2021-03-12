#!/usr/bin/python

import sys
#sys.path.append('/marconi/home/userexternal/cmalosso/thermocepstrum/themocepstrum')
import thermocepstrum as tc
import numpy as np
import argparse
from scipy import signal

def autocorr(x):
    result = signal.correlate(x, x, mode='full', method='fft')
    v = [result[i]/( len(x)-abs( i - (len(x))  ) ) for i in range(len(result))]
    return v[int(result.size/2):]

def green_kubo(data, tmax, dt_fs):
    jen = data[:, 3:]
    temp = data[:, 2]

    Nstep = int(np.rint(tmax / (dt_fs * 1e-3)))
    maxrows = np.size(jen, 0)
    Ncurrs = maxrows // Nstep

    result = np.zeros((Ncurrs, 2))
    for ij in range(Ncurrs):
        init = Nstep * ij
        end = Nstep * (ij + 1) if Nstep * (ij + 1) < jen.shape[0] else jen.shape[0]

        tmean = np.mean(temp[init:end])

        corr = []
        for i in range(3):
            corr.append(np.array(autocorr(jen[init:end, i])))
        mean_corr = (corr[0]+corr[1]+corr[2])/3
        np.savetxt('{}_auto.out'.format(ij), mean_corr)
        result[ij] = [tmean, np.sum(mean_corr[:np.size(mean_corr)//4])]

    return result

def block_analysis(data, tmax, dt_fs, vol):
    jen = data[:, 3:]
    temp = data[:, 2]

    Nstep = int(np.rint(tmax / (dt_fs * 1e-3)))
    maxrows = np.size(jen, 0)
    Ncurrs = maxrows // Nstep

    psds1 = 0
    psds2 = 0

    result = np.zeros([Ncurrs,2])

    for ij in range(Ncurrs):
        init = Nstep * ij
        end = Nstep * (ij + 1) if Nstep * (ij + 1) < jen.shape[0] else jen.shape[0]

        Tmean = np.mean(temp[init:end])

        jj = tc.HeatCurrent(j=jen[init:end], DT_FS=dt_fs, TEMPERATURE=Tmean, units='metal_vis', VOLUME=vol,
                            PSD_FILTER_W=0.05)

        result[ij] = [Tmean, jj.psd[0]*jj.kappa_scale/2]

        # S(w=0) average
        psds1 += (jj.psd[0] - psds1) / (ij + 1)
        psds2 += (jj.psd[0] ** 2 - psds2) / (ij + 1)
    if Ncurrs > 1:
        psdstd = np.sqrt((psds2 - psds1 ** 2) / (Ncurrs - 1))
    else:
        psdstd = 0
    return result, psds1 * jj.kappa_scale / 2, psdstd * jj.kappa_scale / 2

def main():
    parser = argparse.ArgumentParser(description='Block analysis script')
    parser.add_argument('--dt',
                        type=float,
                        required=True,
                        help='time step in femtoseconds')
    parser.add_argument('-b','--box',
                        type=float,
                        required=True,
                        help='box unit size')
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
    vol = args.box**3

    lista = []
    for i in range(first,last+1):
        if args.dt % 1 == 0:
            lista.append(np.loadtxt('stress{}.{:.0f}fs.out'.format(i,args.dt),skiprows=1))
        else:
            lista.append(np.loadtxt('stress{}.{:.1f}fs.out'.format(i, args.dt), skiprows=1))
    data = np.concatenate(lista)


    for time in args.timelist:
        block_result, mean, mean_std = block_analysis(data, time, args.dt, vol)
        print(mean, mean_std)
        np.savetxt('{}ps_result.out'.format(time), block_result)

    kB = 1.3806504
    bartoPa = 1
    for time in args.timelist:

        gk_result = green_kubo(data, time, args.dt)
        for i in range(np.size(gk_result,0)):
            conv = vol / kB / gk_result[i,0] * bartoPa ** 2 * 1e-11
            gk_result[i,1] *= conv
        np.savetxt('{}ps_gk-result.out'.format(time), gk_result)


if __name__ == "__main__":
    main()


