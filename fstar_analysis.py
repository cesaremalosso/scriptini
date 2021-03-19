#!/usr/bin/python

import sys
sys.path.append('/marconi/home/userexternal/cmalosso/thermocepstrum/themocepstrum')
import thermocepstrum as tc
import numpy as np
import argparse


def main():
    parser = argparse.ArgumentParser(description='fstar analysis script')
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
    parser.add_argument('-u', '--units',
                        type = str,
                        required=True,
                        help='units')
    parser.add_argument('-aic', '--aic',
                        type = float,
                        nargs='+',
                        required=True,
                        help='list of AIC coefficients')
    parser.add_argument('-fstar', '--fstar',
                        type = float,
                        nargs='+',
                        required=True,
                        help='list of fstar')
    parser.add_argument('--dir',
                        type=str,
                        required=False,
                        default='./',
                        help='directory of the stress time series')
    args = parser.parse_args()

    first = args.first
    last = args.last
    vol = args.box**3
    time_step = args.dt
    aic_factors = args.aic
    fstar_list = args.fstar

    lista = []
    print('start reading the time series...')
    for i in range(first,last+1):
        if args.dt % 1 == 0:
            lista.append(np.loadtxt('{}stress{}.{:.0f}fs.out'.format(args.dir,i,time_step),skiprows=1))
        else:
            lista.append(np.loadtxt('{}stress{}.{:.1f}fs.out'.format(args.dir,i,time_step), skiprows=1))
    data = np.concatenate(lista)

    jen = data[:,3:6]
    temp = data [:,2]

    if args.units == 'metal':
        u = 'metal_vis'
    elif args.units == 'real':
        u = 'real_vis'

    nnjen=tc.HeatCurrent(jen,units=u,DT_FS=time_step,TEMPERATURE=np.mean(temp),VOLUME=vol,PSD_FILTER_W= 0.1)
    jjjen = {}
    for aic in aic_factors:
        TSKIP_LIST = np.array([nnjen.Nyquist_f_THz/j for j in fstar_list], dtype=np.int)
        jjjen[aic] = tc.heatcurrent.fstar_analysis(nnjen, TSKIP_LIST, Kmin_corrfactor=aic, plot=False)

    output='fstar_analysis.out'

    with open(output,'w') as fs:
        fs.write('#fstar(THz) \t viscosity (cP) \t standard deviation (cP) \n')
        kappa_Kmin={}
        kappa_Kmin_err={}
        for aic in aic_factors:
            FSTAR_LIST = np.array([j.Nyquist_f_THz for j in jjjen[aic]])
            kappa_Kmin[aic] = np.array([j.kappa_Kmin * 100 for j in jjjen[aic]])
            kappa_Kmin_err[aic] = np.array([j.kappa_Kmin_std * 100 for j in jjjen[aic]])

        results = np.copy(FSTAR_LIST)
        for aic in aic_factors:
            results = np.vstack((results,kappa_Kmin[aic],kappa_Kmin_err[aic]))

        np.savetxt(fs,np.transpose(results))

if __name__ == "__main__":
    main()
