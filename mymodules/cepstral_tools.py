import numpy as np
import matplotlib.pyplot as plt
import thermocepstrum as tc

def k_fstar(nnjen, interval=range(1, 20, 2), corrfactors=np.arange(1, 3), plot = False):
    jjjen = {}
    for cor in corrfactors:
        TSKIP_LIST = np.array([nnjen.Nyquist_f_THz / j for j in interval], dtype=int)
        jjjen[cor] = tc.heatcurrent.fstar_analysis(nnjen, TSKIP_LIST, Kmin_corrfactor=cor, plot=False)
    if plot == True:
        plot_k(jjjen, nnjen, TSKIP_LIST, corrfactors, 'kappa vs fstar')
    return jjjen, TSKIP_LIST


def plot_k(jjjen, nnjen, TSKIP_LIST, corrs=np.arange(1,3), title=None):
    f, ax = plt.subplots(1, figsize=(8.0, 6.0), constrained_layout=True)
    ls = 12
    #
    ax.tick_params(axis='x', labelsize=ls)
    ax.tick_params(axis='y', labelsize=ls)
    kappa_Kmin = {}
    kappa_Kmin_err = {}
    Pstar_Kmin = {}
    FSTAR_LIST = nnjen.Nyquist_f_THz / TSKIP_LIST

    for cor in corrs:
        kappa_Kmin[cor] = np.array([j.kappa_Kmin for j in jjjen[cor]])
        kappa_Kmin_err[cor] = np.array([j.kappa_Kmin_std for j in jjjen[cor]])
        Pstar_Kmin[cor] = np.array([j.dct.aic_Kmin + 1 for j in jjjen[cor]])
        f1 = ax.plot(FSTAR_LIST, kappa_Kmin[cor], '-o', label='c={}'.format(cor))
        ax.fill_between(FSTAR_LIST, kappa_Kmin[cor] - kappa_Kmin_err[cor], kappa_Kmin[cor] + kappa_Kmin_err[cor],
                        alpha=0.4)
    if title is not None:
        ax.set_title(title, fontsize=ls)
    ax.set_xlabel('F* (THz)', fontsize=ls)
    ax.set_ylabel(r'$\kappa$ (W/m/K)', fontsize=ls)
    ax.legend(loc='lower right', fontsize=ls, ncol=len(corrs))
    return


def block_analysis(jen, temp, tmax, dt, vol, fstar, corrs=np.arange(1, 3), u='metal_vis'):
    mean = np.zeros(np.size(corrs))
    std = np.zeros(np.size(corrs))
    mean_std = np.zeros(np.size(corrs))

    i = 0
    vis = {}
    vis_std = {}
    for cor in corrs:

        Nstep = int(np.rint(tmax / (dt * 1e-3)))
        maxrows = np.size(jen, 0)
        Ncurrs = maxrows // Nstep
        vis[cor] = []
        vis_std[cor] = []
        t = []
        for ij in range(Ncurrs):
            init = Nstep * ij
            end = Nstep * (ij + 1) if Nstep * (ij + 1) < jen.shape[0] else jen.shape[0]

            tmean = np.mean(temp[init:end])
            t.append(tmean)

            jj = tc.HeatCurrent(j=jen[init:end], DT_FS=dt, TEMPERATURE=tmean, units=u, VOLUME=vol, PSD_FILTER_W=0.3)
            rj = jj.resample_current(fstar_THz=fstar, plot=False, PSD_FILTER_W=0.10)
            rj.cepstral_analysis(Kmin_corrfactor=cor)

            vis[cor].append(rj.kappa_Kmin * 100)
            vis_std[cor].append(rj.kappa_Kmin_std * 100)
        mean[i] = np.mean(vis[cor])
        std[i] = np.std(vis[cor])
        mean_std[i] = np.mean(vis_std[cor])
        i += 1
    return vis, vis_std