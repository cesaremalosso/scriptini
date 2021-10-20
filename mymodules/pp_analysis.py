import numpy as np
import matplotlib.pyplot as plt

a_b=0.529177249
c=29979245800

def read_gdr(filein,output,nmol,vol,dr):
    gdr = np.loadtxt(filein,skiprows=7)
    x = gdr[:,1] * dr
    c = []
    c.append( 8 * np.pi * nmol[0] * nmol[1] / 3 / vol * ((x + dr) ** 3 - x ** 3) )
    c.append( 4 * np.pi * nmol[0] ** 2 / 3 / vol * ((x + dr) ** 3 - x ** 3) )
    c.append( 4 * np.pi * nmol[1] ** 2 / 3 / vol * ((x + dr) ** 3 - x ** 3) )
    gdr_n = np.transpose( np.vstack( (x, gdr[:,6]/c[1], np.sqrt(gdr[:,7])/c[1], gdr[:,4]/c[2], np.sqrt(gdr[:,5])/c[2], gdr[:,2]/c[0], np.sqrt(gdr[:,3])/c[0]) ) )
    with open(output,'w') as fo:
        np.savetxt(output, gdr_n, fmt = '%.18e', header = 'r [Ang] \t g_11 \t sigma_11 \t g_22  \t  sigma_22 \t g_12 \t sigma_12')
    return gdr_n

def plot_gdr(gdr):
    f,ax = plt.subplots(1, figsize=(8.0,6.0), constrained_layout=True)
    ls=12
    #
    ax.tick_params(axis='x', labelsize=ls)
    ax.tick_params(axis='y', labelsize=ls)
    col = ['red', 'blue', 'green']
    lab = ['g_11', 'g_22', 'g_12']
    for i in range(3):
        ax.fill_between(gdr[:,0],gdr[:,2*i+1]-gdr[:,2*i+2],gdr[:,2*i+1]+gdr[:,2*i+2], color = col[i])
        ax.plot(gdr[:,0],gdr[:,2*i+1], label = lab[i], color = col[i])
    plt.legend()

def read_vdos(filein,output,dt,nstep):
    vdos = np.loadtxt(filein)
    x = vdos[:,0] / nstep / dt
    vdos_n = np.zeros((vdos.shape[0],5))
    conv = dt / nstep / 2 / 3
    vdos_n[:,0] = x
    vdos_n[:,1] = conv * (vdos[:,1] + vdos[:,3] + vdos[:,5])
    vdos_n[:,2] = np.sqrt( conv ** 2 * (vdos[:,2] + vdos[:,4] + vdos[:,6]) )
    vdos_n[:,3] = conv * (vdos[:,7] + vdos[:,9] + vdos[:,11])
    vdos_n[:,4] = np.sqrt( conv ** 2 * (vdos[:,8] + vdos[:,10] + vdos[:,12]) )

    with open(output, 'w') as fo:
        np.savetxt(output,vdos_n, fmt = '%.18e', header = 'freq [THz] \t VDOS_1 \t sigma_1 \t VDOS_2 \t  sigma_2')

    return vdos_n

def read_msd(filein,output,dt):
    msd = np.loadtxt(filein, skiprows=7)
    time = [i*dt for i in range(np.shape(msd)[0])]
    msd_n = np.transpose( np.vstack( (time, msd[:,0], np.sqrt(msd[:,1]), msd[:,2], np.sqrt(msd[:,3]),
                                      msd[:,4], np.sqrt(msd[:,5]), msd[:,6], np.sqrt(msd[:,7])) ) )
    with open(output, 'w') as fo:
        np.savetxt(output, msd_n, fmt = '%.18e', header = 'time [ps] \t msd_1 \t sigma_1 \t msd_2 \t sigma_2 \t msd_cm1'
                                                          ' \t sigma_cm1 \t sigma_cm2 \t msd_cm2 ')
    return msd_n

def plot_msd(msd):
    f,ax = plt.subplots(1, figsize=(8.0,6.0), constrained_layout=True)
    ls=12
    #
    ax.tick_params(axis='x', labelsize=ls)
    ax.tick_params(axis='y', labelsize=ls)
    col = ['red', 'blue']
    lab = ['msd_1', 'msd_2']
    for i in range(2):
        ax.fill_between(msd[:,0],msd[:,2*i+1]-msd[:,2*i+2],msd[:,2*i+1]+msd[:,2*i+2], color = col[i])
        ax.plot(msd[:,0],msd[:,2*i+1], label = lab[i], color = col[i])
    plt.legend()
