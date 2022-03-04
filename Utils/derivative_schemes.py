"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Derivative schemes functions

@description: Derivative schemes used for post-processing.
            Here, the python function np.gradient is used and corrected at the edges.
            Those schemes are NOT using the non-physical point at the end of each direction, i.e.
            with MULTIFAST language, all arrays have a shape of (n3-1)*(n2-1)*(n1-1)
"""

import numpy as np

def D1_3Dy(f, Yc, f_wall = 0.0, correct_at_wall = True):
    """Compute the 1st order derivative of a 3D or 1D array in the vertical direction.
    
    Parameters
    ----------
    f : 3D numpy array with shape (nz,ny,nx) OR 1D numpy array with shape (ny)
        Fields to be derivated
    Yc : 1D numpy array with shape (ny)
        Position of the nodes in the vertical direction (cell-centered)
    f_wall : Real
        Value at the wall of the input field 
    correct_at_wall : Boolean
        Trigger the option to correct the 1st order derivative at the wall by using the value at the wall
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx) OR 1D numpy array with shape (ny)
        The vertical 1st order derivative of the filed given as input
            
    """
    
    # The flag correct_at_wall is used only if values at wall are available
    
    if (len(f.shape)>2):
		# That means we compute the derivative of a 3D array

	    # dimensions of array
        n3  = f.shape[0]
        n2  = f.shape[1]
        n1  = f.shape[2]
        
        d1f = np.zeros([n3,n2,n1])
        a_down = np.zeros([3])
        
        d1f = np.gradient(f, Yc, edge_order=2, axis=1)
        
        if (correct_at_wall):
        
            al_down=(Yc[1]-Yc[0])/(Yc[0])
            
            a_down[0]=-al_down/(al_down+1)
            a_down[1]=(al_down-1)/al_down
            a_down[2]=1/(al_down**2+al_down)
    
            a_down=a_down/( (Yc[0]) )
    
            d1f[:,0,:]=(a_down[0]*0.0 + a_down[1]*f[:,0,:]+ a_down[2]*f[:,1,:])
        
    else:
		# That means we compute the derivative of a 1D array

	    # dimensions of array
        n2  = f.shape[0]
        
        d1f = np.zeros([n2])
        a_down = np.zeros([3])
        
        d1f = np.gradient(f, Yc, edge_order=2)
        
        if (correct_at_wall):
        
            al_down=(Yc[1]-Yc[0])/(Yc[0])
            
            a_down[0]=-al_down/(al_down+1)
            a_down[1]=(al_down-1)/al_down
            a_down[2]=1/(al_down**2+al_down)
    
            a_down=a_down/( (Yc[0]) )
    
            d1f[0]=(a_down[0]*0.0 + a_down[1]*f[0]+ a_down[2]*f[1])
        
    return(d1f)



def D2_3Dy(f, Yc, Y):
    """Compute the 2nd order derivative of a 3D or 1D array in the vertical direction.
    
    Parameters
    ----------
    f : 3D numpy array with shape (nz,ny,nx) OR 1D numpy array with shape (ny)
        Fields to be derivated
    Yc : 1D numpy array with shape (ny)
        Position of the nodes in the vertical direction (cell-centered)
    Y : 1D numpy array with shape (ny)
        Position of the nodes in the vertical direction (staggered)
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx) OR 1D numpy array with shape (ny)
        The vertical 2nd order derivative of the filed given as input
            
    """
    
    if (len(f.shape)>2):
        # That means we compute the derivative of a 3D array
        
        # dimensions of array
        n3  = f.shape[0]
        n2  = f.shape[1]
        n1  = f.shape[2]
        
        d2f = np.zeros([n3,n2,n1])
        
        d2f = np.gradient( np.gradient(f, Yc, edge_order=2, axis=1), \
                       Yc, edge_order=2, axis=1)
            
        dx2 = 2/n2
        beta = ( (Yc[0]-Y[0])**2 ) / dx2**2
        alpha = (Yc[1]-Yc[0]) / (Yc[0]-Y[0])
        a1_d = 2 / ((alpha+1)*beta)
        a2_d = -2 / (alpha*beta)
        a3_d = 2 / ((alpha**2+alpha)*beta)
        
        d2f[:,0,:] = ( a1_d*0 + a2_d*f[:,0,:] + a3_d*f[:,1,:] ) / dx2**2
        
    else:
		# That means we compute the derivative of a 1D array

	    # dimensions of array
        n2  = f.shape[0]
        
        d2f = np.zeros([n2])
        
        d2f = np.gradient( np.gradient(f, Yc, edge_order=2), \
                       Yc, edge_order=2)
            
        dx2 = 2/n2
        beta = ( (Yc[0]-Y[0])**2 ) / dx2**2
        alpha = (Yc[1]-Yc[0]) / (Yc[0]-Y[0])
        a1_d = 2 / ((alpha+1)*beta)
        a2_d = -2 / (alpha*beta)
        a3_d = 2 / ((alpha**2+alpha)*beta)
        
        d2f[0] = ( a1_d*0 + a2_d*f[0] + a3_d*f[1] ) / dx2**2
        
    return(d2f)