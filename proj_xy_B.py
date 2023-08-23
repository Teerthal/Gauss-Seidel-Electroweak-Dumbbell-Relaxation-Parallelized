##Python scripts to read and plot magnetic field data in the xy plane##
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
ini_list = ['%s/%s' % (master_path,file)
            for file in os.listdir(master_path) if file.startswith('B_xy_')
            and file.endswith('dat')]

arr = []

for file in ini_list:
    data = np.genfromtxt(file)
    [arr.append(i) for i in data]

arr = np.array(arr)

itimes = np.unique(arr[:, 0]).astype(int)
xuniques = np.unique(arr[:, 1]).astype(int)
yuniques = np.unique(arr[:, 2]).astype(int)
z_slices = np.unique(arr[:, 3]).astype(int)


def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"

def vector_field_plot():
    for z_idx, z_slc in enumerate(z_slices):
        fig,axs = plt.subplots(ncols=2)
        for idx,i in enumerate([itimes[0],itimes[-2]]):
            ax=axs[idx]
            stride=6

            tar = np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                                np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                                ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0], tar)

            X = arr[tar, 1]
            Y = arr[tar, 2]
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]
            B_z = -arr[tar, 6]
            mod_B = np.sqrt(B_x**2+B_y**2)
            B_x = -arr[tar, 4]/mod_B
            B_y = -arr[tar, 5]/mod_B

            df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
            xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
            df = df.sort_values(by=[xcol, ycol])
            xvals = df[xcol].unique()
            yvals = df[ycol].unique()
            Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
            By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

            #plt.quiver(xvals,yvals,Bx_vals,By_vals)
            ax.quiver(X,Y,B_x,B_y,scale=40.0)
            ax.set_title('i:%i'%i)
            #ax.set_box_aspect(1)
            if idx==0:
                ax.set_xlabel('x')
                ax.set_ylabel('y')
        fig.savefig('%s/%s_xy_vector_%i_z_%i.pdf'%(plot_path,'B',i,z_slc))
    return
vector_field_plot()

def vector_scaled_plot():
    for z_idx, z_slc in enumerate(z_slices):
        fig,axs = plt.subplots(ncols=2)
        for idx,i in enumerate([itimes[0],itimes[-2]]):
            ax=axs[idx]
            stride=4

            tar = np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                                np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                                ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0], tar)

            X = arr[tar, 1]
            Y = arr[tar, 2]
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]
            B_z = -arr[tar, 6]
            mod_B = np.sqrt(B_x**2+B_y**2)
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]

            df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
            xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
            df = df.sort_values(by=[xcol, ycol])
            xvals = df[xcol].unique()
            yvals = df[ycol].unique()
            Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
            By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

            #plt.quiver(xvals,yvals,Bx_vals,By_vals)
            ax.quiver(X,Y,B_x,B_y,scale_units='xy',scale=0.01,width=0.0035,headwidth=2.5,pivot='middle')
            ax.set_title('i:%i'%i)
            #ax.set_box_aspect(1)
            ax.set_xlim([-50,50])
            ax.set_ylim([-50,50])
            if idx==0:
                ax.set_xlabel('x')
                ax.set_ylabel('y')
        fig.savefig('%s/%s_xy_vector_unscaled_%i_z_%i.pdf'%(plot_path,'B',i,z_slc))
    return
vector_scaled_plot()


def stream_field_plot():
    for z_idx, z_slc in enumerate(z_slices):
        fig,axs = plt.subplots(ncols=2)
        for idx,i in enumerate([itimes[0],itimes[-2]]):
            ax=axs[idx]

            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0],
             np.argwhere(arr[:, 0] == i)[:, 0])

            X = arr[tar, 1]
            Y = arr[tar, 2]
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]
            B_z = -arr[tar, 6]
            mod_B = np.sqrt(B_x**2+B_y**2)
            B_x = -arr[tar, 4]/mod_B
            B_y = -arr[tar, 5]/mod_B

            df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
            xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
            df = df.sort_values(by=[xcol, ycol])
            xvals = df[xcol].unique()
            yvals = df[ycol].unique()
            Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
            By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T
            ax.streamplot(xvals,yvals,Bx_vals,By_vals,density=2,linewidth=0.5,
            arrowsize=0.5,arrowstyle='->')
            #ax.streamplot(X,Y,B_x,B_y)
            ax.set_title('i:%i'%i)
            #ax.set_box_aspect(1)
            if idx==0:
                ax.set_xlabel('x')
                ax.set_ylabel('y')
        fig.savefig('%s/%s_xy_stream_%i_z_%i.pdf'%(plot_path,'B',i,z_slc))
    return
stream_field_plot()


def mag_contours():
    
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        
        for z_idx, z_slc in enumerate(z_slices):
            fig,axs = plt.subplots()
            ax=axs
            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0]
                       , np.argwhere(arr[:, 0] == i)[:, 0])
                       
            X = arr[tar, 1]
            Y = arr[tar, 2]
            Z = np.sqrt(arr[tar, 4]**2+arr[tar, 5]**2+arr[tar, 6]**2)
            #lvls=np.logspace(1e-2,10,50)s
            #im=plt.tricontourf(X, Y, Z,levels=lvls)#,locator=ticker.LogLocator())
            #cntrs=plt.tricontour(X, Y, Z,[0.01,0.1,0.5,1.0,1.5,2.0],colors='white',linewidths=1.0)

            lev_exp = np.linspace(np.floor(np.log10(Z.min())),
                            np.ceil(np.log10(Z.max())+1),50)
            levs = np.power(10, lev_exp)
            #im = plt.tricontourf(X, Y, Z, levs, norm=colors.LogNorm())
            im = ax.tricontourf(X, Y, Z, levs,
                                norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))

            #plt.clabel(cntrs,cntrs.levels,inline=True,fontsize=10,fmt='%1.1f')
            
            #ax.set_title('z:%i'%z_slc)
            #plt.xlim([-75,75])
            #plt.ylim([-75,75])
            #plt.colorbar(im)
            ##ax.set_box_aspect(1)
            cbar = plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
            cbar.ax.set_yticks(levs[::5])
            cbar.ax.set_yticklabels(["{:.2e}".format(tck) for tck in levs[::5]])
            # This is the fix for the white lines between contour levels
            #if idx==0:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            for c in im.collections:
                c.set_edgecolor("face")
    
            plt.tight_layout()
            fig.savefig('%s/%s_xy_mag_contours_%i_z_%i.pdf'%(plot_path,'B',i,z_slc))
    return
mag_contours()

def mag_cntr_stream():
    
    for idx,i in enumerate([itimes[0],itimes[-2]]):
        
        for z_idx, z_slc in enumerate(z_slices):
            fig,axs = plt.subplots()
            ax=axs
            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0]
                       , np.argwhere(arr[:, 0] == i)[:, 0])

            dx = 0.1                       
            X = arr[tar, 1]*dx
            Y = arr[tar, 2]*dx
            Z = np.sqrt(arr[tar, 4]**2+arr[tar, 5]**2+arr[tar, 6]**2)
            #lvls=np.logspace(1e-2,10,50)s
            #im=plt.tricontourf(X, Y, Z,levels=lvls)#,locator=ticker.LogLocator())
            #cntrs=plt.tricontour(X, Y, Z,[0.01,0.1,0.5,1.0,1.5,2.0],colors='white',linewidths=1.0)

            lev_exp = np.linspace(np.floor(np.log10(Z.min())),
                            np.ceil(np.log10(Z.max())+1),75)
            levs = np.power(10, lev_exp)
            #im = plt.tricontourf(X, Y, Z, levs, norm=colors.LogNorm())
            im = ax.tricontourf(X, Y, Z, levs,
                                norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max())
                                ,cmap='plasma')
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]
            B_z = -arr[tar, 6]
            mod_B = np.sqrt(B_x**2+B_y**2)
            B_x = -arr[tar, 4]/mod_B
            B_y = -arr[tar, 5]/mod_B
            df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
            xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
            df = df.sort_values(by=[xcol, ycol])
            xvals = df[xcol].unique()
            yvals = df[ycol].unique()
            Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
            By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T
            ax.streamplot(xvals,yvals,Bx_vals,By_vals,density=2,linewidth=0.65,
            arrowsize=0.75,arrowstyle='->',color='k')
            #plt.clabel(cntrs,cntrs.levels,inline=True,fontsize=10,fmt='%1.1f')
            
            #ax.set_title('z:%i'%z_slc)
            #plt.xlim([-75,75])
            #plt.ylim([-75,75])
            #plt.colorbar(im)
            ##ax.set_box_aspect(1)
            cbar = plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
            cbar.ax.set_yticks(levs[::5])
            cbar.ax.set_yticklabels(["{:.1e}".format(tck) for tck in levs[::5]])
            # This is the fix for the white lines between contour levels
            #if idx==0:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            for c in im.collections:
                c.set_edgecolor("face")
            ax.set_xlim([-150*dx,150*dx])
            ax.set_ylim([-150*dx,150*dx])
            plt.tight_layout()
            fig.savefig('%s/%s_xy_mag_cntrs_stream_%i_z_%i.pdf'%(plot_path,'B',i,z_slc))
    return
mag_cntr_stream()

def combined_plot():
    for z_idx, z_slc in enumerate(z_slices):
        for idx,i in enumerate([itimes[0],itimes[-2]]):
            fig = plt.figure()
            ax=fig.add_subplot()


            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0]
            , np.argwhere(arr[:, 0] == i)[:, 0])
                       
            X = arr[tar, 1]
            Y = arr[tar, 2]
            Z = np.sqrt(arr[tar, 4]**2+arr[tar, 5]**2+arr[tar, 6]**2)
            lev_exp = np.linspace(np.floor(np.log10(Z.min())),
                            np.ceil(np.log10(Z.max())+1),50)
            levs = np.power(10, lev_exp)
            im = ax.tricontourf(X, Y, Z, levs,
                                norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()))
            cbar = plt.colorbar(im, ax=ax)#,fraction=0.046, pad=0.04)
            cbar.ax.set_yticks(levs[::5])
            cbar.ax.set_yticklabels(["{:.1e}".format(tck) for tck in levs[::5]])
            for c in im.collections:
                c.set_edgecolor("face")


            stride=4

            tar = np.intersect1d(np.intersect1d(np.argwhere(arr[:, 0] == i)[:, 0],
                                np.argwhere(np.isin(arr[:, 1],xuniques[::stride]))[:, 0])
                                ,np.argwhere(np.isin(arr[:, 2],yuniques[::stride]))[:, 0])
            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0], tar)

            X = arr[tar, 1]
            Y = arr[tar, 2]
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]
            B_z = -arr[tar, 6]
            mod_B = np.sqrt(B_x**2+B_y**2)
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]

            df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y))
            xcol, ycol, bxcol,bycol = "x", "y", "bx", "by"
            df = df.sort_values(by=[xcol, ycol])
            xvals = df[xcol].unique()
            yvals = df[ycol].unique()
            Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
            By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T

            #plt.quiver(xvals,yvals,Bx_vals,By_vals)
            ax.quiver(X,Y,B_x,B_y,scale_units='xy',scale=0.03,width=0.0035,headwidth=2.5,pivot='middle')
            #ax.set_title('i:%i'%i)
            #ax.set_box_aspect(1)
            ax.set_xlim([-50,50])
            ax.set_ylim([-50,50])
            #if idx==0:
            ax.set_xlabel('x')
            ax.set_ylabel('y')

            fig.savefig('%s/%s_xy_combined_%i_z_%i.pdf'%(plot_path,'B',i,z_slc))

    return

combined_plot()



def stream_unscaled_plot():
    for z_idx, z_slc in enumerate(z_slices):
        fig,axs = plt.subplots(ncols=2)
        for idx,i in enumerate([itimes[0],itimes[-2]]):
            ax=axs[idx]

            tar = np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0],
             np.argwhere(arr[:, 0] == i)[:, 0])
            
            xy_mask=20
            # tar = np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(arr[:, 3] == z_slc)[:, 0],
            # np.argwhere(arr[:, 0] == i)[:, 0]), np.argwhere(np.abs(arr[:, 1]) > xy_mask)[:, 0])
            # ,np.argwhere(np.abs(arr[:, 2]) > xy_mask)[:, 0])
            arr[np.sqrt(np.abs(arr[:,1])**2+np.abs(arr[:,2])**2)<xy_mask,4]=0
            arr[np.sqrt(np.abs(arr[:,1])**2+np.abs(arr[:,2])**2)<xy_mask,5]=0

            X = arr[tar, 1]
            Y = arr[tar, 2]
            B_x = -arr[tar, 4]
            B_y = -arr[tar, 5]
            B_z = -arr[tar, 6]
            mod_B = np.sqrt(B_x**2+B_y**2)
            B_x = -arr[tar, 4]#/mod_B
            B_y = -arr[tar, 5]#/mod_B

            df = pd.DataFrame(dict(x=X, y=Y, bx=B_x,by=B_y,bz=B_z))
            xcol, ycol, bxcol,bycol,bzcol = "x", "y", "bx", "by", "bz"
            df = df.sort_values(by=[xcol, ycol])
            xvals = df[xcol].unique()
            yvals = df[ycol].unique()
            Bx_vals = df[bxcol].values.reshape(len(xvals), len(yvals)).T
            By_vals = df[bycol].values.reshape(len(xvals), len(yvals)).T
            Bz_vals = df[bzcol].values.reshape(len(xvals), len(yvals)).T
            lw=5*np.sqrt(Bx_vals**2+By_vals**2+Bz_vals**2)/np.max(np.sqrt(Bx_vals**2+By_vals**2+Bz_vals**2))
            ax.streamplot(xvals,yvals,Bx_vals,By_vals,density=1,linewidth=0.5,
            arrowsize=0.5,arrowstyle='->')
            #ax.streamplot(X,Y,B_x,B_y)
            #ax.set_title('i:%i'%i)
            #ax.set_box_aspect(1)
            
            ax.set_xlabel('x')
            ax.set_ylabel('y')
        fig.savefig('%s/%s_xy_stream_unscaled_%i_z_%i.pdf'%(plot_path,'B',i,z_slc))
    return
stream_unscaled_plot()
