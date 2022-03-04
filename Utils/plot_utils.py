"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Utilities functions for 2D plotting

@description: Contains general functions for 2D plotting in the three configurations (XY,XZ,ZY).
"""

import numpy as np

from CFD_plot import CFD_plot
from CFD_mesh import CFD_mesh
from Settings import Settings


from utilities import find_nearest, print_infos


def plot_2D_XY(arrayToShow, zValue, mesh:CFD_mesh, settings:Settings):
    """Plot 2D fields in XY configuration using the CFD_plot class.
    
    Take a slice of the 3D array 'arrayToShow' at the position 'zValue'.
    The index is obtained with the fonction find_nearest(mesh.Zc, zValue)
    
    Parameters
    ----------
    arrayToShow : Numpy array of shape (nz,ny,nx)
        Array to plot
    zValue : Real
        Position of the slice to plot
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    settings : Settings
        Settings of the simulation
            
    """

    fig = CFD_plot('full')

    if settings.contours:
        fig.add_contour(mesh.ZY,mesh.YZ,arrayToShow[find_nearest(mesh.Zc, zValue),:,:], cmap='seismic')

    if settings.colored_contours:
        fig.add_contourf(mesh.ZY,mesh.YZ,arrayToShow[find_nearest(mesh.Zc, zValue),:,:], cmap='seismic')

    fig.chg_x_axis(r'$X$',axis_low_bound=0,axis_high_bound=mesh.Zc[mesh.nz-1]+(mesh.Zc[1]-mesh.Zc[0])/2)
    fig.chg_y_axis(r'$Y$',axis_low_bound=0,axis_high_bound=2,axis_ticks=[0,1,2])
    fig.custom_layout()
    fig.display()

def plot_2D_XZ(arrayToShow, yValue, mesh:CFD_mesh, settings:Settings):
    """Plot 2D fields in XZ configuration using the CFD_plot class.
    
    Take a slice of the 3D array 'arrayToShow' at the position 'yValue'.
    The index is obtained with the fonction find_nearest(mesh.Yc, yValue)
    
    Parameters
    ----------
    arrayToShow : Numpy array of shape (nz,ny,nx)
        Array to plot
    zValue : Real
        Position of the slice to plot
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    settings : Settings
        Settings of the simulation
            
    """

    fig = CFD_plot('full')

    if settings.contours:
        fig.add_contour(mesh.ZX,mesh.XZ,np.transpose(arrayToShow[:,find_nearest(mesh.Yc, yValue),:]), cmap='seismic')

    if settings.colored_contours:
        fig.add_contourf(mesh.ZX,mesh.XZ,np.transpose(arrayToShow[:,find_nearest(mesh.Yc, yValue),:]), cmap='seismic')

    fig.chg_x_axis(r'$X$',axis_low_bound=0,axis_high_bound=mesh.Zc[mesh.nz-1]+(mesh.Zc[1]-mesh.Zc[0])/2)
    fig.chg_y_axis(r'$Z$',axis_low_bound=mesh.Xc[0]-mesh.X[mesh.nx//2 + 1],axis_high_bound=mesh.Xc[-1]-mesh.X[mesh.nx//2 + 1])
    fig.custom_layout()
    fig.display()

def plot_2D_ZY(arrayToShow, xValue, mesh:CFD_mesh, settings:Settings):
    """Plot 2D fields in ZY configuration using the CFD_plot class.
    
    Take a slice of the 3D array 'arrayToShow' at the position 'xValue'.
    The index is obtained with the fonction find_nearest(mesh.Xc, xValue)
    
    Parameters
    ----------
    arrayToShow : Numpy array of shape (nz,ny,nx)
        Array to plot
    zValue : Real
        Position of the slice to plot
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    settings : Settings
        Settings of the simulation
            
    """

    fig = CFD_plot('full')

    if settings.contours:
        fig.add_contour(mesh.ZY,mesh.YZ,arrayToShow[:,:,find_nearest(mesh.Xc, xValue)], cmap='seismic')

    if settings.colored_contours:
        fig.add_contourf(mesh.ZY,mesh.YZ,arrayToShow[:,:,find_nearest(mesh.Xc, xValue)], cmap='seismic')

    fig.chg_x_axis(r'$Z$',axis_low_bound=mesh.Xc[0]-mesh.X[mesh.nx//2 + 1],axis_high_bound=mesh.Xc[-1]-mesh.X[mesh.nx//2 + 1])
    fig.chg_y_axis(r'$Y$',axis_low_bound=0,axis_high_bound=2,axis_ticks=[0,1,2])
    fig.custom_layout()
    fig.display()

# def plot_2D_XY_TE(arrayToShow, xValue, axe):

#     if contours:
#         fig.add_contour(ZY,YZ,arrayToShow[find_nearest(Xc, zValue),:,:], cmap='binary')

#     if colored_contours:
#         fig.add_contourf(ZY,YZ,arrayToShow[find_nearest(Xc, zValue),:,:], cmap='seismic')

# def plot_2D_ZY_TE(arrayToShow, zValue, axe):
#     tmp = max(np.max(arrayToShow[find_nearest(Xc, zValue),:,:]), abs(np.min(arrayToShow[find_nearest(Xc, zValue),:,:])))
#     colorbar = np.linspace(-tmp,tmp,40)

#     if contours:
#         axe.contour(ZY,YZ,arrayToShow[find_nearest(Xc, zValue),:,:], colors='black')
#         axe.set_xlabel('Z',fontsize=lsize)

#     if colored_contours:
#         axe.contourf(ZY,YZ,arrayToShow[find_nearest(Xc, zValue),:,:], colorbar, cmap='seismic')
#         axe.set_xlabel('Z',fontsize=lsize)

# def plot_2D_XZ_TE(arrayToShow, yValue, axe):
#     tmp = max(np.max(arrayToShow[:,find_nearest(Yc, yValue),:]), abs(np.min(arrayToShow[:,find_nearest(Yc, yValue),:])))
#     colorbar = np.linspace(-tmp,tmp,40)
#     # colorbar = np.arange(-tmp,tmp,1.5*10**(-6))

#     if contours:
#         axe.contour(ZX,XZ,np.transpose(arrayToShow[:,find_nearest(Yc, yValue),:]), colorbar, colors='black')
#         axe.set_xlabel('X',fontsize=lsize)

#     if colored_contours:
#         axe.contourf(ZX,XZ,np.transpose(arrayToShow[:,find_nearest(Yc, yValue),:]), colorbar, cmap='seismic')
#         axe.set_xlabel('X',fontsize=lsize)  
    
def plot_2D_contours(array3D, arrayName, mesh:CFD_mesh, settings:Settings):
    """Plot 2D fields in any configuration using the CFD_plot class.
    
    Depending settings of user, call function plotting 2D field in ZY, XY, or XZ configuration
    
    Parameters
    ----------
    array3D : Numpy array of shape (nz,ny,nx)
        Name of the file to write
    arrayName : String
        Name of the array to plot
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    settings : Settings
        Settings of the simulation
            
    """

    if (type(settings.plot_x_equal_to) != bool):

        plot_2D_ZY(array3D, settings.plot_x_equal_to)
        print_infos(array3D[:,:,find_nearest(mesh.Zc, settings.plot_x_equal_to)], arrayName)

    if (type(settings.plot_y_equal_to) != bool):

        plot_2D_XZ(array3D, settings.plot_y_equal_to)
        print_infos(array3D[:,find_nearest(mesh.Yc, settings.plot_y_equal_to),:], arrayName)

    if (type(settings.plot_z_equal_to) != bool):

        plot_2D_XY(array3D, settings.plot_z_equal_to)
        print_infos(array3D[find_nearest(mesh.Xc, settings.plot_z_equal_to),:,:], arrayName)
        

def plot_2D_contours_with_spot(array2D, spot2D, mesh:CFD_mesh, settings:Settings):
    """Plot 2D fields with the spot in XZ configuration using the CFD_plot class.
    
    Parameters
    ----------
    array2D : Numpy array of shape (nz,nx)
        Array to plot
    spot2D : Numpy array of shape (nz,nx)
        Spot to plot
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    settings : Settings
        Settings of the simulation
            
    """

    fig = CFD_plot('full')
    fig.add_contourf(mesh.XZ,mesh.ZX,array2D, cmap='seismic')
    fig.add_contour(mesh.XZ,mesh.ZX,spot2D, cmap='binary', linewidths=2)
    fig.chg_x_axis(r'$X$',axis_low_bound=mesh.Xc[0],axis_high_bound=mesh.Xc[-1])
    fig.chg_y_axis(r'$Z$',axis_low_bound=0-settings.Lz/2,axis_high_bound=settings.Lz-settings.Lz/2)
    fig.custom_layout()
    fig.display()