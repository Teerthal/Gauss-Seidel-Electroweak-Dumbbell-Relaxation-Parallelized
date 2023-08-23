import numpy as np
import matplotlib.pyplot as plt
import os as os
from numpy import reshape
import plotly.graph_objects as go
from matplotlib import ticker,cm,colors
import matplotlib.tri as tri
import sys
master_path = '/path/to/directory'

plot_path = master_path+'/images'

ini_list = ['%s/%s' %(master_path, file)
            for file in os.listdir(master_path) if file.startswith('phi_xz_')
            and file.endswith('dat')]


arr = []

for file in ini_list:
    data = np.genfromtxt(file)
    [arr.append(i) for i in data]

arr = np.array(arr)
itimes = np.unique(arr[:, 0]).astype(int)

dx=0.1


def phi_contours_xz():

    for idx,i in enumerate([itimes[0],itimes[-2]]):
        fig,axs = plt.subplots()
        ax=axs
        tar = np.argwhere(arr[:, 0] == i)[:, 0]
        
        X = arr[tar, 1]
        Y = arr[tar, 2]
        Z = arr[tar, 3]

        #im = plt.tricontourf(X, Y, Z, levs, norm=colors.LogNorm())
        lev_exp = np.linspace(np.log10(Z.min()),np.log10(Z.max()),75)
        levs = np.power(10, lev_exp)
        #plt.clabel(cntrs,cntrs.levels,inline=True,fontsize=10,fmt='%1.1f')
        # im = ax.tricontourf(X, Y, Z, levs,
        #                 norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max())
        #                 ,cmap='plasma')
        im = ax.tricontourf(X*dx, Y*dx, Z,cmap='plasma_r',levels=20,vmax=0.95)
        #ax.set_title('z:%i'%z_slc)
        plt.xlim([-100*dx,100*dx])
        plt.ylim([-100*dx,100*dx])
        cbar=plt.colorbar(im)
        cbar.set_ticks([0.00,0.25,0.50,0.75,0.95])
        ax.set_box_aspect(1)
        # cbar = plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
        # cbar.ax.set_yticks(levs[::5])
        # cbar.ax.set_yticklabels(["{:.2e}".format(tck) for tck in levs[::5]])
        # This is the fix for the white lines between contour levels
        #if idx==0:
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        for c in im.collections:
            c.set_edgecolor("face")

        fig.tight_layout()
        fig.savefig('%s/xz_phi_contours_%i.pdf'%(plot_path,i))
        plt.close()
    return
phi_contours_xz()
