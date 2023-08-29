"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Interpolation schemes functions

@description: Interpolation schemes based on the DRP numerical scheme of MULTIFAST.
            More details can be found in the paper of Bauer (2015, doi: 10.1016/j.compfluid.2014.10.009)
            Here, the scheme are interpolating the data from the face center (where the data is in the output of MULTIFAST) 
            to the cell center (for post-processing and visualization).
            Those schemes are NOT using the non-physical point at the end of each direction, i.e.
            with MULTIFAST language, all arrays have a shape of (n3-1)*(n2-1)*(n1-1)
"""

import numpy as np

def D0s_DRP5_3Dz(f):
    """Interpolate the given 3D field in the spanwise direction from a staggered configuration to a cell-centered one.
    
    Parameters
    ----------
    f : 3D numpy array with shape (nz,ny,nx)
        Fields to be interpolated
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The interpolated 3D field, with the variable at the cell-center
            
    """
    
    # dimensions of array
    n3  = f.shape[0]
    n2  = f.shape[1]
    n1  = f.shape[2]

    ff = np.zeros([n3,n2,n1])
                    
    for i in range(n3-1):
        ff[i,:,:] = (f[i+1,:,:]+f[i,:,:])/2
        
    ff[n3-1,:,:] = (f[0,:,:]+f[n3-1,:,:])/2
    
                            
    return(ff)


def D0s_DRP5_3Dy(f, f_wall = 0.0):
    """Interpolate the given 3D field in the vertical direction from a staggered configuration to a cell-centered one.
    
    Parameters
    ----------
    f : 3D numpy array with shape (nz,ny,nx)
        Fields to be interpolated
    f_wall : Real
        Value at the wall of the field
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The interpolated 3D field, with the variable at the cell-center
            
    """
    
    # f_wall is assumed to be 0 in default case, as for most of the case in MULTIFAST
    
    # dimensions of array
    n3  = f.shape[0]
    n2  = f.shape[1]
    n1  = f.shape[2]

    ff = np.zeros([n3,n2,n1])
                                                
    for i in range(n2-1):
        ff[:,i,:] = (f[:,i+1,:]+f[:,i,:])/2
        
    ff[:,n2-1,:] = (f_wall+f[:,n2-1,:])/2

    return(ff)


def D0s_DRP5_3Dx(f):
    """Interpolate the given 3D field in the streamwise direction from a staggered configuration to a cell-centered one.
    
    Parameters
    ----------
    f : 3D numpy array with shape (nz,ny,nx)
        Fields to be interpolated
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The interpolated 3D field, with the variable at the cell-center
            
    """
    
    # dimensions of array
    n3  = f.shape[0]
    n2  = f.shape[1]
    n1  = f.shape[2]

    ff = np.zeros([n3,n2,n1])
                                    
    for i in range(n1-1):
        ff[:,:,i] = (f[:,:,i+1]+f[:,:,i])/2
        
    ff[:,:,n1-1] = (f[:,:,0]+f[:,:,n1-1])/2

    return(ff)
