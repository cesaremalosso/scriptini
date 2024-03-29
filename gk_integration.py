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

def green_kubo(jen, temp, tmax, gk_tmax, dt_fs):

    Nstep = int(np.rint(tmax / (dt_fs * 1e-3)))
    maxrows = np.size(jen, 0)
    Ncurrs = maxrows // Nstep

    result = np.zeros((Ncurrs, 8))
    for ij in range(Ncurrs):
        init = Nstep * ij
        end = Nstep * (ij + 1) if Nstep * (ij + 1) < jen.shape[0] else jen.shape[0]

        tmean = np.mean(temp[init:end])*191/192

        corr = []
        for i in range(3):
            corr.append(np.array(autocorr(jen[init:end, i])))
        mean_corr = (corr[0]+corr[1]+corr[2])/3
        np.savetxt('{}_auto.out'.format(ij), mean_corr)
        result[ij] = [tmean, np.sum(mean_corr[:int(25e3)]), np.sum(mean_corr[:int(50e3)]), np.sum(mean_corr[:int(75e3)]),
                      np.sum(mean_corr[:int(100e3)]), np.sum(mean_corr[:int(200e3)]), np.sum(mean_corr[:int(300e3)]), np.sum(mean_corr[:int(400e3)])]

    return result

def block_analysis(jen, temp, tmax, dt_fs, vol):

    Nstep = int(np.rint(tmax / (dt_fs * 1e-3)))
    maxrows = np.size(jen, 0)
    Ncurrs = maxrows // Nstep

    psds1 = 0
    psds2 = 0

    result = np.zeros([Ncurrs,3])

    for ij in range(Ncurrs):
        init = Nstep * ij
        end = Nstep * (ij + 1) if Nstep * (ij + 1) < jen.shape[0] else jen.shape[0]

        Tmean = np.mean(temp[init:end])*191/192

        jj = tc.HeatCurrent(j=jen[init:end], DT_FS=dt_fs, TEMPERATURE=Tmean, units='metal_vis', VOLUME=vol,
                            PSD_FILTER_W=0.05)

        rj = jj.resample_current(fstar_THz=4, plot=False, PSD_FILTER_W=0.10)
        rj.cepstral_analysis(Kmin_corrfactor=2)

        result[ij] = [Tmean, rj.dct.psd[0]*jj.kappa_scale/2, rj.kappa_Kmin_std]

        # S(w=0) average
        psds1 += (rj.dct.psd[0] - psds1) / (ij + 1)
        psds2 += (rj.dct.psd[0] ** 2 - psds2) / (ij + 1)
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
    parser.add_argument('--tmax',
                        type = int,
                        required=True,
                        help='green-kubo extremum of integration in ps')
    parser.add_argument('--timelist',
                        type = int,
                        nargs='+',
                        required=True,
                        help='list of the size of blocks for block-analysis')
    parser.add_argument('-u', '--units',
                        type = str,
                        required=True,
                        help='units real or metal supported')
    args = parser.parse_args()

    first = args.first
    last = args.last
    vol = args.box**3
    units = args.units

    lista = []
    for i in range(first,last+1):
        if args.dt % 1 == 0:
            lista.append(np.loadtxt('stress{}.{:.0f}fs.out'.format(i,args.dt),skiprows=1))
        else:
            lista.append(np.loadtxt('stress{}.{:.1f}fs.out'.format(i, args.dt), skiprows=1))
    data = np.concatenate(lista)

    jen = data[:,3:6]
    temp = data [:,2]

    if units == 'metal':
        print('metal units')
    elif units == 'real':
        jen *= 1.01325
        print('real units')

    for time in [100,200,500]:
        block_result, mean, mean_std = block_analysis(jen, temp, time, args.dt, vol)
        print(mean, mean_std)
        with open('{}ps_result.out'.format(time),'w') as f:
            np.savetxt(f, block_result)
            f.write('\t \t {} \n'.format(np.mean(block_result[:,2])))
            f.write('# mean values \n')
            f.write('{} \t {} \t {}'.format(np.mean(block_result[:,0]),mean,mean_std))

    kB = 1.3806504
    bartoPa = 1
    for time in args.timelist:
        gk_result = green_kubo(jen, temp, time, args.tmax, args.dt)
        for i in range(np.size(gk_result,0)):
            conv = vol / kB / gk_result[i,0] * args.dt * bartoPa ** 2 * 1e-11
            gk_result[i,1:] *= conv
        with open('{}ps_gk-result.out'.format(time),'w') as f:
#            f.write('# mean temperature \t tmax = {} \t tmax = {} \t tmax = {}\n'.format(args.tmax/2,
#                                                                                       args.tmax, args.tmax*2))
            np.savetxt(f, gk_result)
            mean = np.array([np.mean(gk_result[:,i]) for i in range(8)])
            f.write('# mean values \n')
            np.savetxt(f, mean.reshape(1, mean.shape[0]))
#            f.write('{} \t {} \t {} \t {}'.format(np.mean(gk_result[:,0]), np.mean(gk_result[:,1]), np.mean(gk_result[:,2]), np.mean(gk_result[:,3])))


if __name__ == "__main__":
    main()


