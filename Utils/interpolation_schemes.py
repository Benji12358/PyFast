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
    
    A = np.zeros(6)
    
    A[1] =  1.234102595369109     /2
    A[2] = -0.3184044667712       /2
    A[3] =  0.11029870162898      /2
    A[4] = -0.030619166038291     /2
    A[5] =  0.0046223358114003    /2
    

    for i in range(4,n3-5):
        ff[i, :, :]= A[5]*(f[i+5, :, :] + f[i-4, :, :]) \
                   + A[4]*(f[i+4, :, :] + f[i-3, :, :]) \
                   + A[3]*(f[i+3, :, :] + f[i-2, :, :]) \
                   + A[2]*(f[i+2, :, :] + f[i-1, :, :]) \
                   + A[1]*(f[i+1, :, :] + f[i-0, :, :])

    # Assuming the z direction is periodic !
    for i in range(4):
        ff[i, :, :]=  A[5]*(f[i+5, :, :] + f[(n3+i-4)%(n3), :, :])  \
                    + A[4]*(f[i+4, :, :] + f[(n3+i-3)%(n3), :, :])  \
                    + A[3]*(f[i+3, :, :] + f[(n3+i-2)%(n3), :, :])  \
                    + A[2]*(f[i+2, :, :] + f[(n3+i-1)%(n3), :, :])  \
                    + A[1]*(f[i+1, :, :] + f[(n3+i-0)%(n3), :, :]) 

    for i in range(n3-5,n3):
        ff[i, :, :]=  A[5]*(f[(i+5)%(n3), :, :] + f[i-4, :, :]) \
                    + A[4]*(f[(i+4)%(n3), :, :] + f[i-3, :, :]) \
                    + A[3]*(f[(i+3)%(n3), :, :] + f[i-2, :, :]) \
                    + A[2]*(f[(i+2)%(n3), :, :] + f[i-1, :, :]) \
                    + A[1]*(f[(i+1)%(n3), :, :] + f[i-0, :, :])
                            
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
    
    A = np.zeros(6)
    
    A[1] =  1.234102595369109    /2
    A[2] = -0.3184044667712      /2
    A[3] =  0.11029870162898     /2
    A[4] = -0.030619166038291    /2
    A[5] =  0.0046223358114003   /2

    for j in range(4,n2-5):
        ff[:, j, :]=  A[5]*(f[:, j+5, :] + f[:, j-4, :]) \
                    + A[4]*(f[:, j+4, :] + f[:, j-3, :]) \
                    + A[3]*(f[:, j+3, :] + f[:, j-2, :]) \
                    + A[2]*(f[:, j+2, :] + f[:, j-1, :]) \
                    + A[1]*(f[:, j+1, :] + f[:, j, :])

    # Assuming the y direction has walls !
    # A dirichlet condition is then used
    ff[:,0,:] = 0.5*(f[:,1,:]+f[:,0,:])
            
    ff[:,1,:]  = 9.0/16*(f[:,2,:]+f[:,1,:]) \
                -0.0625*(f[:,3,:]+f[:,0,:])
                    
    ff[:,2,:]  = 0.59395104312381* (f[:,3,:]+f[:,2,:]) \
                -0.10967656468571* (f[:,4,:]+f[:,1,:]) \
                +0.015725521561903*(f[:,5,:]+f[:,0,:])
            
    ff[:,3,:]  = 0.59395104312381* (f[:,4,:]+f[:,3,:]) \
                -0.10967656468571* (f[:,5,:]+f[:,2,:]) \
                +0.015725521561903*(f[:,6,:]+f[:,1,:])
            
    ff[:,n2-1,:] = 0.5*(f[:,n2-1,:] + f_wall)        

    ff[:,n2-2,:] = 9.0/16*(f[:,n2-1,:]+f[:,n2-2,:]) \
                -0.0625*(f_wall+f[:,n2-3,:])
            
    ff[:,n2-3,:] = 0.59395104312381 * (f[:,n2-2,:] + f[:,n2-3,:]) \
                -  0.10967656468571 * (f[:,n2-1,:] + f[:,n2-4,:]) \
                +  0.015725521561903* (f_wall + f[:,n2-5,:])
                            
    ff[:,n2-4,:] = 0.59395104312381 * (f[:,n2-3,:] + f[:,n2-4,:]) \
                -  0.10967656468571 * (f[:,n2-2,:] + f[:,n2-5,:]) \
                +  0.015725521561903* (f[:,n2-1,:] + f[:,n2-6,:])
                
    ff[:,n2-5,:] = A[1] * (f[:,n2-4,:] + f[:,n2-5,:]) \
                +  A[2] * (f[:,n2-3,:] + f[:,n2-6,:]) \
                +  A[3] * (f[:,n2-2,:] + f[:,n2-7,:]) \
                +  A[4] * (f[:,n2-1,:] + f[:,n2-8,:]) \
                +  A[5] * (f_wall      + f[:,n2-9,:])

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
    
    A = np.zeros(6)

    A[1] =  1.234102595369109    /2
    A[2] = -0.3184044667712      /2
    A[3] =  0.11029870162898     /2
    A[4] = -0.030619166038291    /2
    A[5] =  0.0046223358114003   /2


    for k in range(4,n1-5):
        ff[:, :, k]=  A[5]*(f[:, :, k+5] + f[:, :, k-4]) \
                    + A[4]*(f[:, :, k+4] + f[:, :, k-3]) \
                    + A[3]*(f[:, :, k+3] + f[:, :, k-2]) \
                    + A[2]*(f[:, :, k+2] + f[:, :, k-1]) \
                    + A[1]*(f[:, :, k+1] + f[:, :, k-0])

    # Assuming the x direction is periodic !
    for k in range(4):
        ff[:, :, k]=  A[5]*(f[:, :, k+5] + f[:, :, (n1+k-4)%(n1)]) \
                    + A[4]*(f[:, :, k+4] + f[:, :, (n1+k-3)%(n1)]) \
                    + A[3]*(f[:, :, k+3] + f[:, :, (n1+k-2)%(n1)]) \
                    + A[2]*(f[:, :, k+2] + f[:, :, (n1+k-1)%(n1)]) \
                    + A[1]*(f[:, :, k+1] + f[:, :, (n1+k-0)%(n1)])

    for k in range(n1-5,n1):
        ff[:, :, k]=  A[5]*(f[:, :, (k+5)%(n1)] + f[:, :, k-4]) \
                    + A[4]*(f[:, :, (k+4)%(n1)] + f[:, :, k-3]) \
                    + A[3]*(f[:, :, (k+3)%(n1)] + f[:, :, k-2]) \
                    + A[2]*(f[:, :, (k+2)%(n1)] + f[:, :, k-1]) \
                    + A[1]*(f[:, :, (k+1)%(n1)] + f[:, :, k-0])

    return(ff)