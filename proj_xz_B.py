##Python scripts to read and plot magnetic field data in the xz plane##
#Has various magnetic field plotting functions that can be used as desired#

import numpy as np
import matplotlib.pyplot as plt
import os as os
import sys
from numpy import reshape
import matplotlib.tri as tri
from matplotlib import ticker,cm,colors
# import dash
# import dash_core_components as dcc
# import dash_html_components as html
import pandas as pd

master_path='/path/to/directory'

plot_path = master_path+'/images'

#Collects all data files with the correct designations#
#Due to the parallelized nature, each parallel
#domain dumps data into different files marked by process number

ini_list = ['%s/%s' % (master_path,file)
            for file in os.listdir(master_path) if file.startswith('B_xz_')
            and file.endswith('dat')]

arr = []

for file in ini_list:
    data = np.genfromtxt(file)
    [arr.append(i) for i in data]

arr = np.array(arr)
print(np.shape(arr))
itimes = np.unique(arr[:, 0]).astype(int)
xuniques = np.unique(arr[:, 1]).astype(int)
yuniques = np.unique(arr[:, 2]).astype(int)
print('itimes:%s'%itimes)
print(xuniques)
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"


def vector_field_plot():
    fig,axs = plt.subplots(ncols=2)
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        ax=axs[idx]
        stride=6
        tar = np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                            np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                            ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
        X = arr[tar, 1]
        Y = arr[tar, 2]
        B_x = -arr[tar, 3]
        B_y = -arr[tar, 4]
        B_z = -arr[tar, 5]
        mod_B = np.sqrt(B_x**2+B_y**2+B_z**2)
        B_x = -arr[tar, 3]/mod_B
        B_y = -arr[tar, 4]/mod_B
        B_z = -arr[tar, 5]/mod_B

        df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
        xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
        df = df.sort_values(by=[xcol, ycol])
        xvals = df[xcol].unique()
        yvals = df[ycol].unique()
        Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
        By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

        #plt.quiver(xvals,yvals,Bx_vals,By_vals)
        ax.quiver(X,Y,B_x,B_z,scale=40.0)
        ax.set_title('i:%i'%i)
        # ax.set_box_aspect(1)
        if idx==0:
            ax.set_xlabel('x')
            ax.set_ylabel('z')
    fig.savefig('%s/%s_xz_vector_%i.pdf'%(plot_path,'B',i))
    plt.close()
    return
vector_field_plot()

def vector_unscaled_plot():
    #fig,axs = plt.subplots(ncols=2)
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        fig = plt.figure()
        ax=fig.add_subplot()
        #ax=axs[idx]
        stride=6
        tar = np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                            np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                            ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
        X = arr[tar, 1]
        Y = arr[tar, 2]
        B_x = -arr[tar, 3]
        B_y = -arr[tar, 4]
        B_z = -arr[tar, 5]
        mod_B = np.sqrt(B_x**2+B_y**2+B_z**2)
        B_x = -arr[tar, 3]
        B_y = -arr[tar, 4]
        B_z = -arr[tar, 5]

        df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
        xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
        df = df.sort_values(by=[xcol, ycol])
        xvals = df[xcol].unique()
        yvals = df[ycol].unique()
        Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
        By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

        #plt.quiver(xvals,yvals,Bx_vals,By_vals)
        ax.quiver(X,Y,B_x,B_z,scale_units='xy',scale=0.01,width=0.0035,headwidth=2.5,pivot='middle')
        ax.set_title('i:%i'%i)
        # ax.set_box_aspect(1)
        if idx==0:
            ax.set_xlabel('x')
            ax.set_ylabel('z')
        ax.set_xlim([-150,150])
        ax.set_ylim([-150,150])
        fig.savefig('%s/%s_xz_vector_unscaled_%i.pdf'%(plot_path,'B',i))
        plt.close()
    return
vector_unscaled_plot()

def stream_field_plot():
    fig,axs = plt.subplots(ncols=2)
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        ax=axs[idx]
        tar = np.argwhere(arr[:, 0] == i)[:, 0]
        X = arr[tar, 1]
        
        Y = arr[tar, 2]
        B_x = -arr[tar, 3]
        B_y = -arr[tar, 4]
        B_z = -arr[tar, 5]
        mod_B = np.sqrt(B_x**2+B_y**2+B_z**2)
        B_x = -arr[tar, 3]/mod_B
        B_y = -arr[tar, 4]/mod_B
        B_z = -arr[tar, 5]/mod_B

        df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_z))
        xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
        df = df.sort_values(by=[xcol, ycol])
        xvals = df[xcol].unique()
        yvals = df[ycol].unique()
        Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
        By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

        #plt.quiver(xvals,yvals,Bx_vals,By_vals)
        #ax.streamplot(X,Y,B_x,B_z)
        ax.streamplot(xvals,yvals,Bx_vals,By_vals,density=3,linewidth=0.5,
        arrowsize=0.5,arrowstyle='->')
        #ax.set_title('i:%i'%i)
        # ax.set_box_aspect(1)
        if idx==0:
            ax.set_xlabel('x')
            ax.set_ylabel('z')
    fig.savefig('%s/%s_xz_stream_%i.pdf'%(plot_path,'B',i))
    plt.close()
    return
stream_field_plot()

def mag_contours():
    fig,axs = plt.subplots(ncols=2)
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        ax=axs[idx]
        tar = np.argwhere(arr[:, 0] == i)[:, 0]
        X = arr[tar, 1]
        Y = arr[tar, 2]
        Z = np.sqrt(arr[tar, 3]**2+arr[tar, 4]**2+arr[tar, 5]**2)
        #lvls=np.logspace(1e-2,10,50)
        #im=plt.tricontourf(X, Y, Z,levels=lvls)#,locator=ticker.LogLocator())
        #cntrs=plt.tricontour(X, Y, Z,[0.01,0.1,0.5,1.0,1.5,2.0],colors='white',linewidths=1.0)


        #im = plt.tricontourf(X, Y, Z, levs, norm=colors.LogNorm())
        
        ###########Log scale contours###############
        lev_exp = np.linspace(np.floor(np.log10(Z.min())),
                           np.ceil(np.log10(Z.max())+1),50)
        levs = np.power(10, lev_exp)

        im = ax.tricontourf(X, Y, Z, levs,
                            norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
        ###########Log scale contours###############

        #plt.clabel(cntrs,cntrs.levels,inline=True,fontsize=10,fmt='%1.1f')
        
        ax.set_title('i:%i'%i)
        #plt.xlim([-75,75])
        #plt.ylim([-75,75])
        #plt.colorbar(im)
        # ax.set_box_aspect(1)
        cbar = plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
        cbar.ax.set_yticks(levs[::5])
        cbar.ax.set_yticklabels(["{:.2e}".format(tck) for tck in levs[::5]])
        # This is the fix for the white lines between contour levels
        if idx==0:
            ax.set_xlabel('x')
            ax.set_ylabel('z')
        for c in im.collections:
            c.set_edgecolor("face")
     
        
    plt.tight_layout()
    fig.savefig('%s/%s_xz_mag_contours_%i.pdf'%(plot_path,'B',i))
    plt.close()
    return
mag_contours()


def combined_plot():
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        fig = plt.figure()
        ax=fig.add_subplot()

        tar = np.argwhere(arr[:, 0] == i)[:, 0]
        dx = 0.1
        X = arr[tar, 1]*dx
        Y = arr[tar, 2]*dx
        Z = np.sqrt(arr[tar, 3]**2+arr[tar, 4]**2+arr[tar, 5]**2)
        ###########Log scale contours###############
        lev_exp = np.linspace(np.floor(np.log10(Z.min())),
                           np.ceil(np.log10(Z.max())+1),75)
        levs = np.power(10, lev_exp)

        im = ax.tricontourf(X, Y, Z, levs, cmap='plasma',
                            norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
        ###########Log scale contours###############        
        # ax.set_box_aspect(1)
        cbar = plt.colorbar(im, ax=ax)#,fraction=0.046, pad=0.04)
        cbar.ax.set_yticks(levs[::5])
        cbar.ax.set_yticklabels(["{:.1e}".format(tck) for tck in levs[::5]])
        for c in im.collections:
            c.set_edgecolor("face")
     
        tar = np.argwhere(arr[:, 0] == i)[:, 0]
        
        dx = 0.1
        X = arr[tar, 1]*dx
        
        Y = arr[tar, 2]*dx
        # B_x = -arr[tar, 3]
        # B_y = -arr[tar, 4]
        # B_z = -arr[tar, 5]

        B_x = -arr[tar, 3]
        B_y = -arr[tar, 4]
        B_z = -arr[tar, 5]
        mod_B = np.sqrt(B_x**2+B_y**2+B_z**2)        
        # modB_mask=1e-1
        # B_x[mod_B<modB_mask]=0
        # B_y[mod_B<modB_mask]=0
        # B_z[mod_B<modB_mask]=0

        B_x = B_x/mod_B
        B_y = B_y/mod_B
        B_z = B_z/mod_B

        df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_z,bz=B_y,mb=mod_B))
        xcol, ycol, bxcol,bzcol,bycol,mbcol = "x", "y", "bx", "by", "bz","mb"
        df = df.sort_values(by=[xcol, ycol])
        xvals = df[xcol].unique()
        yvals = df[ycol].unique()
        Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
        By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T
        Bz_vals = df[bzcol].values.reshape(len(xvals), len(yvals)).T
        mB_vals = df[mbcol].values.reshape(len(xvals), len(yvals)).T
        lw=1.5*(np.log10(mB_vals)*1.2+np.abs(np.min(np.log10(mB_vals)))
        )/np.abs(np.min(np.log10(mB_vals)))
        clrs=np.log10(mB_vals)
        ax.streamplot(xvals,yvals,Bx_vals,Bz_vals,density=2.5,linewidth=0.65,
        arrowsize=0.5,arrowstyle='->',color='k',maxlength=400)
        # ax.streamplot(xvals,yvals,Bx_vals,Bz_vals,density=1.6,linewidth=lw,
        # arrowsize=0.5,arrowstyle='->',color=clrs,cmap='gist_yarg')
        ax.set_xlim([-150*dx,150*dx])
        ax.set_ylim([-150*dx,150*dx])
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        fig.savefig('%s/%s_xz_combined_%i.pdf'%(plot_path,'B',i))
    return
combined_plot()

def combined_plot2():
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        fig = plt.figure()
        ax=fig.add_subplot()

        tar = np.argwhere(arr[:, 0] == i)[:, 0]
        X = arr[tar, 1]
        Y = arr[tar, 2]
        Z = np.sqrt(arr[tar, 3]**2+arr[tar, 4]**2+arr[tar, 5]**2)
        ###########Log scale contours###############
        lev_exp = np.linspace(np.floor(np.log10(Z.min())),
                           np.ceil(np.log10(Z.max())+1),50)
        levs = np.power(10, lev_exp)

        im = ax.tricontourf(X, Y, Z, levs,
                            norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
        ###########Log scale contours###############        
        # ax.set_box_aspect(1)
        cbar = plt.colorbar(im, ax=ax)#,fraction=0.046, pad=0.04)
        cbar.ax.set_yticks(levs[::5])
        cbar.ax.set_yticklabels(["{:.1e}".format(tck) for tck in levs[::5]])
        for c in im.collections:
            c.set_edgecolor("face")

        stride=4
        tar = np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                            np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                            ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
        X = arr[tar, 1]
        Y = arr[tar, 2]
        B_x = -arr[tar, 3]
        B_y = -arr[tar, 4]
        B_z = -arr[tar, 5]
        mod_B = np.sqrt(B_x**2+B_y**2+B_z**2)
        B_x = -arr[tar, 3]/mod_B
        B_y = -arr[tar, 4]/mod_B
        B_z = -arr[tar, 5]/mod_B

        df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
        xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
        df = df.sort_values(by=[xcol, ycol])
        xvals = df[xcol].unique()
        yvals = df[ycol].unique()
        Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
        By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

        #plt.quiver(xvals,yvals,Bx_vals,By_vals)
        ax.quiver(X,Y,B_x,B_z,scale=40.0,width=0.002,headwidth=1,headaxislength=4,headlength=7)

        ax.set_xlabel('x')
        ax.set_ylabel('z')
        fig.savefig('%s/%s_xz_combined2_%i.pdf'%(plot_path,'B',i))
    return
combined_plot2()


def combined_plot2_x():
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        fig = plt.figure()
        ax=fig.add_subplot()

        # tar = np.argwhere(arr[:, 0] == i)[:, 0]
        # X = arr[tar, 1]
        # Y = arr[tar, 2]
        # Z = np.sqrt(arr[tar, 3]**2+arr[tar, 4]**2+arr[tar, 5]**2)
        # ###########Log scale contours###############
        # lev_exp = np.linspace(np.floor(np.log10(Z.min())),
        #                    np.ceil(np.log10(Z.max())+1),50)
        # levs = np.power(10, lev_exp)

        # im = ax.tricontourf(X, Y, Z, levs,
        #                     norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
        # ###########Log scale contours###############        
        # ax.set_box_aspect(1)
        # cbar = plt.colorbar(im, ax=ax)#,fraction=0.046, pad=0.04)
        # cbar.ax.set_yticks(levs[::5])
        # cbar.ax.set_yticklabels(["{:.1e}".format(tck) for tck in levs[::5]])
        # for c in im.collections:
        #     c.set_edgecolor("face")

        stride=4
        tar = np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                            np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                            ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
                            ,np.argwhere(arr[:, 2] == 0)[:, 0])
        X = arr[tar, 1]
        Y = arr[tar, 2]

        # B_x = -arr[tar, 3]
        # B_y = -arr[tar, 4]
        # B_z = -arr[tar, 5]
        # mod_B = np.sqrt(B_x**2+B_y**2+B_z**2)
        # B_x = -arr[tar, 3]/mod_B
        # B_y = -arr[tar, 4]/mod_B
        # B_z = -arr[tar, 5]/mod_B

        B_x = arr[tar, 3]
        B_y = arr[tar, 4]
        B_z = arr[tar, 5]
        # mod_B = np.sqrt(B_x**2+B_y**2+B_z**2)
        # B_x = arr[tar, 3]/mod_B
        # B_y = arr[tar, 4]/mod_B
        # B_z = arr[tar, 5]/mod_B


        df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
        xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
        df = df.sort_values(by=[xcol, ycol])
        xvals = df[xcol].unique()
        yvals = df[ycol].unique()
        Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
        By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

        #plt.quiver(xvals,yvals,Bx_vals,By_vals)
        ax.quiver(X,Y,B_x,B_z,scale=5.0,width=0.002,headwidth=4,headaxislength=4,headlength=4)

        ax.set_xlabel('x')
        ax.set_ylabel('z')
        fig.savefig('%s/%s_xz_combined2_x_%i.pdf'%(plot_path,'B',i))
    return
combined_plot2_x()

def B_x_dependence():
    for idx,i in enumerate(itimes):
        fig = plt.figure()
        ax=fig.add_subplot()

        stride=1
        tar = np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                            np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                            ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
                            ,np.argwhere(arr[:, 2] == 0)[:, 0])
        
        
        X = arr[tar, 1]
        B_x = arr[tar, 3]
        B_y = arr[tar, 4]
        B_z = arr[tar, 5]
        B_mag = np.sqrt(B_x**2+B_y**2+B_z**2)
        ax.scatter(X,B_mag,label='|B|',facecolor='')
        ax.scatter(X,np.abs(B_x),label='|B_x|',facecolor='')
        ax.scatter(X,np.abs(B_y),label='|B_y|',facecolor='')
        ax.scatter(X,np.abs(B_z),label='|B_z|',facecolor='')

        sub_arr = arr[tar,:]
        sub_arr.sort(axis=1)

        X = sub_arr[1]
        B_x = sub_arr[3]
        B_y = sub_arr[4]
        B_z = sub_arr[5]

        B_mag = np.sqrt(B_x**2+B_y**2+B_z**2)
        print(np.shape(B_x))
        ax.plot(X,B_mag,label='|B|')
        ax.plot(X,np.abs(B_x),label='|B_x|')
        ax.plot(X,np.abs(B_y),label='|B_y|')
        ax.plot(X,np.abs(B_z),label='|B_z|')
        ax.set_yscale('Log')
        plt.legend()
        fig.savefig('%s/%s_B_x_dependence_%i.pdf'%(plot_path,'B',i))
    return

# B_x_dependence()
def B_z_dependence():
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        fig = plt.figure()
        ax=fig.add_subplot()

        stride=1
        tar = np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                            np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                            ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
                            ,np.argwhere(arr[:, 1] == 0)[:, 0])
        X = arr[tar, 1]
        Y = arr[tar, 2]

        B_x = arr[tar, 3]
        B_y = arr[tar, 4]
        B_z = arr[tar, 5]

        B_mag = np.sqrt(B_x**2+B_y**2+B_z**2)

        ax.scatter(Y,B_mag)
        ax.set_yscale('Log')
        fig.savefig('%s/%s_B_z_dependence_%i.pdf'%(plot_path,'B',i))
    return

# B_z_dependence()

def energy_over_iterations():
    file = master_path +'/total_energy.dat'
    data = np.genfromtxt(file)
    print(np.shape(data))
    fig = plt.figure()
    ax=fig.add_subplot()
    ax.plot(data[:,0],data[:,1])
    ax.set_yscale('Log')
    fig.savefig('%s/%s_Energy_v_iterations.pdf'%(plot_path,'B'))
    return
energy_over_iterations()
