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
from Settings import Settings

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
    
            d1f[:,0,:]=(a_down[0]*f_wall + a_down[1]*f[:,0,:]+ a_down[2]*f[:,1,:])
        
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
    
            d1f[0]=(a_down[0]*f_wall + a_down[1]*f[0]+ a_down[2]*f[1])
        
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

def D1_3Dy_IBM(d1f, f, settings:Settings, Yc, f_wall = 0.0):
    """Compute the 1st order derivative of a 3D array in the vertical direction.
    
    Parameters
    ----------
    d1f : 3D numpy array with shape (nz,ny,nx)
        Fields to be derivated
    Yc : 1D numpy array with shape (ny)
        Position of the nodes in the vertical direction (cell-centered)
    f_wall : Real
        Value at the wall of the input field 
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The vertical 1st order derivative of the field given as input
            
    """
	# That means we compute the derivative of a 3D array
    
    for n in range(len(settings.i_start)):
        
        j_s  = settings.j_start[n]
        j_e  = settings.j_end[n]
        
        if (j_s==0):
            # case roughness on bottom wall

            vert_pos = j_e + 1
            vert_pos_2 = j_e + 2
            
            d1f[:,vert_pos,:]=(f[:,vert_pos_2,:]-f[:,vert_pos,:])/(Yc[vert_pos_2]-Yc[vert_pos])
            
            
            vert_pos = j_e
            vert_pos_2 = j_e - 1
                        
            d1f[:,vert_pos,:]=(f[:,vert_pos,:]-f[:,vert_pos_2,:])/(Yc[vert_pos]-Yc[vert_pos_2]) # gives 0 whatever we set in the object
            
            
            
        else:
            # case roughness on top wall

            vert_pos = j_s - 1
            vert_pos_2 = j_s - 2
            
            d1f[:,vert_pos,:]=(f[:,vert_pos_2,:]-f[:,vert_pos,:])/(Yc[vert_pos_2]-Yc[vert_pos])
            
            
            vert_pos = j_s
            vert_pos_2 = j_s + 1
                        
            d1f[:,vert_pos,:]=(f[:,vert_pos,:]-f[:,vert_pos_2,:])/(Yc[vert_pos]-Yc[vert_pos_2]) # gives 0 whatever we set in the object
                
    return(d1f)

def D2_3Dy_IBM(d2f, f, settings:Settings, Yc, dx2, f_wall = 0.0):
    """Compute the 1st order derivative of a 3D array in the vertical direction.
    
    Parameters
    ----------
    d1f : 3D numpy array with shape (nz,ny,nx)
        Fields to be derivated
    Yc : 1D numpy array with shape (ny)
        Position of the nodes in the vertical direction (cell-centered)
    f_wall : Real
        Value at the wall of the input field 
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The vertical 1st order derivative of the field given as input
            
    """
	# That means we compute the derivative of a 3D array
    
    for n in range(len(settings.i_start)):
        
        i_s  = settings.i_start[n]
        i_e  = settings.i_end[n]
        
        j_s  = settings.j_start[n]
        j_e  = settings.j_end[n]
        
        y_e = settings.y_end[n]
        
        k_s  = settings.k_start[n]
        k_e  = settings.k_end[n]
        
        if (j_s==0):
            
            # pour j_e - 1
            vert_pos = j_e - 1
            vert_pos_aft = j_e
            vert_pos_bef = j_e - 2
            
            hp_over_hm_up=(Yc[vert_pos_aft]-Yc[vert_pos])/(Yc[vert_pos]-Yc[vert_pos_bef])

            hm2_dx2q_u=((Yc[vert_pos]-Yc[vert_pos_bef])**2)/dx2**2

            a1_u    =   2/((hp_over_hm_up+1)*hm2_dx2q_u)
            a2_u    =   -2/(hp_over_hm_up*hm2_dx2q_u)
            a3_u    =   2/((hp_over_hm_up**2+hp_over_hm_up)*hm2_dx2q_u)
            
            d2f[:,vert_pos,:]= (f[:,vert_pos_aft,:]*a3_u+f[:,vert_pos,:]*a2_u+f[:,vert_pos_bef,:]*a1_u)/dx2**2            
            
            # pour j_e
            vert_pos = j_e
            vert_pos_bef = j_e - 1
            vert_pos_bef_bef = j_e - 2
            
            delta_1 = (Yc[vert_pos_bef]-Yc[vert_pos])
            delta_2 = (Yc[vert_pos_bef_bef]-Yc[vert_pos])
            
            d2f[:,vert_pos,:]= 2*(f[:,vert_pos_bef_bef,:]*(-delta_1)+f[:,vert_pos,:]*(delta_1-delta_2)+f[:,vert_pos_bef,:]*(delta_2))/(delta_2*delta_1**2-delta_1*delta_2**2)
            
            # pour j_e+1
            vert_pos = j_e + 1
            vert_pos_aft = j_e + 2
            vert_pos_aft_aft = j_e + 3
            
            delta_1 = (Yc[vert_pos_aft]-Yc[vert_pos])
            delta_2 = (Yc[vert_pos_aft_aft]-Yc[vert_pos])
            
            d2f[:,vert_pos,:]= 2*(f[:,vert_pos_aft_aft,:]*(-delta_1)+f[:,vert_pos,:]*(delta_1-delta_2)+f[:,vert_pos_aft,:]*(delta_2))/(delta_2*delta_1**2-delta_1*delta_2**2)            
            
            # pour j_e + 2
            vert_pos = j_e + 2
            vert_pos_aft = j_e + 1
            vert_pos_bef = j_e + 3
            
            hp_over_hm_down=(Yc[vert_pos_aft]-Yc[vert_pos])/(Yc[vert_pos]-Yc[vert_pos_bef])
            
            hm2_dx2q_d=((Yc[vert_pos]-Yc[vert_pos_bef])**2)/dx2**2
            
            a1_d    =   2/((hp_over_hm_down+1)*hm2_dx2q_d)
            a2_d    =   -2/(hp_over_hm_down*hm2_dx2q_d)
            a3_d    =   2/((hp_over_hm_down**2+hp_over_hm_down)*hm2_dx2q_d)
            
            d2f[:,vert_pos,:]= (f[:,vert_pos_aft,:]*a1_d+f[:,vert_pos,:]*a2_d+f[:,vert_pos_bef,:]*a3_d)/dx2**2
            
        else:
            # vert_pos = j_s - 1
            # shift = -1
            # y_object = settings.y_start[n]
            
            # al_tmp=abs(Yc[vert_pos+shift]-Yc[vert_pos])/abs(Yc[vert_pos]-y_object)
            
            # a_coefs[0]=-al_tmp/(al_tmp+1)
            # a_coefs[1]=(al_tmp-1)/al_tmp
            # a_coefs[2]=1/(al_tmp**2+al_tmp)
    
            # a_coefs=a_coefs/( abs(Yc[vert_pos]-y_object) )
            
            # d1f[k_s:k_e+1,vert_pos,i_s:i_e+1]=(a_coefs[0]*f_wall + a_coefs[1]*f[k_s:k_e+1,vert_pos,i_s:i_e+1]+ a_coefs[2]*f[k_s:k_e+1,vert_pos+shift,i_s:i_e+1])
            a=1
            
    return(d2f)

def D1_3Dx_IBM(d1f, f, settings:Settings, Xc, f_wall = 0.0):
    """Compute the 1st order derivative of a 3D array in the vertical direction.
    
    Parameters
    ----------
    d1f : 3D numpy array with shape (nz,ny,nx)
        Fields to be derivated
    Xc : 1D numpy array with shape (nx)
        Position of the nodes in the streamwise direction (cell-centered)
    f_wall : Real
        Value at the wall of the input field 
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The streamwise 1st order derivative of the field given as input
            
    """
	# That means we compute the derivative of a 3D array
    
    for n in range(len(settings.i_start)):
        
        i_s  = settings.i_start[n]
        i_e  = settings.i_end[n]

        # Downstream face
        vert_pos = i_e + 1
        vert_pos_2 = i_e + 2
        
        d1f[:,:,vert_pos]=(f[:,:,vert_pos_2]-f[:,:,vert_pos])/(Xc[vert_pos_2]-Xc[vert_pos])
        
        
        vert_pos = i_e
        vert_pos_2 = i_e - 1
                    
        d1f[:,:,vert_pos]=(f[:,:,vert_pos]-f[:,:,vert_pos_2])/(Xc[vert_pos]-Xc[vert_pos_2]) # gives 0 whatever we set in the object

        # Upstream face
        vert_pos = i_s - 1
        vert_pos_2 = i_s - 2
        
        d1f[:,:,vert_pos]=(f[:,:,vert_pos_2]-f[:,:,vert_pos])/(Xc[vert_pos_2]-Xc[vert_pos])
        
        
        vert_pos = i_s
        vert_pos_2 = i_s + 1
                    
        d1f[:,:,vert_pos]=(f[:,:,vert_pos]-f[:,:,vert_pos_2])/(Xc[vert_pos]-Xc[vert_pos_2]) # gives 0 whatever we set in the object
                
    return(d1f)

def D1_3Dz_IBM(d1f, f, settings:Settings, Zc, f_wall = 0.0):
    """Compute the 1st order derivative of a 3D array in the vertical direction.
    
    Parameters
    ----------
    d1f : 3D numpy array with shape (nz,ny,nx)
        Fields to be derivated
    Zc : 1D numpy array with shape (nz)
        Position of the nodes in the spanwise direction (cell-centered)
    f_wall : Real
        Value at the wall of the input field 
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The vertical 1st order derivative of the field given as input
            
    """
	# That means we compute the derivative of a 3D array
    
    for n in range(len(settings.i_start)):
        
        k_s  = settings.k_start[n]
        k_e  = settings.k_end[n]
        
        # Downstream face
        vert_pos = (k_e + 1)%len(Zc)
        vert_pos_2 = (k_e + 2)%len(Zc)
        
        d1f[vert_pos,:,:]=(f[vert_pos_2,:,:]-f[vert_pos,:,:])/(Zc[vert_pos_2]-Zc[vert_pos])
        
        
        vert_pos = k_e%len(Zc)
        vert_pos_2 = k_e - 1
                    
        d1f[vert_pos,:,:]=(f[vert_pos,:,:]-f[vert_pos_2,:,:])/(Zc[vert_pos]-Zc[vert_pos_2]) # gives 0 whatever we set in the object 

        # Upstream face
        vert_pos = k_s - 1
        vert_pos_2 = k_s - 2
        
        d1f[vert_pos,:,:]=(f[vert_pos_2,:,:]-f[vert_pos,:,:])/(Zc[vert_pos_2]-Zc[vert_pos])
        
        
        vert_pos = k_s
        vert_pos_2 = k_s + 1
                    
        d1f[vert_pos,:,:]=(f[vert_pos,:,:]-f[vert_pos_2,:,:])/(Zc[vert_pos]-Zc[vert_pos_2]) # gives 0 whatever we set in the object
                
    return(d1f)
        


def D1_1Dy_IBM(d1f, f, settings:Settings, Yc, f_wall = 0.0):
    """Compute the 1st order derivative of a 3D array in the vertical direction.
    
    Parameters
    ----------
    d1f : 3D numpy array with shape (nz,ny,nx)
        Fields to be derivated
    Yc : 1D numpy array with shape (ny)
        Position of the nodes in the vertical direction (cell-centered)
    f_wall : Real
        Value at the wall of the input field 
    
    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The vertical 1st order derivative of the filed given as input
            
    """
	# That means we compute the derivative of a 3D array
    
    a_coefs = np.zeros([3])
    
    for n in range(len(settings.i_start)):
        
        j_s  = settings.j_start[n]
        j_e  = settings.j_end[n]
        
        if (j_s==0):
            vert_pos = j_e + 1
            shift = 1
            y_object = settings.y_end[n] 
            
            al_tmp=abs(Yc[vert_pos]-y_object)/abs(Yc[vert_pos+shift]-Yc[vert_pos])
            
            a_coefs[0]=-al_tmp/(al_tmp+1)
            a_coefs[1]=(al_tmp-1)/al_tmp
            a_coefs[2]=1/(al_tmp**2+al_tmp)
    
            a_coefs=a_coefs/( abs(Yc[vert_pos+shift]-Yc[vert_pos]) )
            
            d1f[vert_pos]=(a_coefs[2]*f_wall + a_coefs[1]*f[vert_pos]+ a_coefs[0]*f[vert_pos+shift])
        else:
            vert_pos = j_s - 1
            shift = -1
            y_object = settings.y_start[n]
            
            al_tmp=abs(Yc[vert_pos+shift]-Yc[vert_pos])/abs(Yc[vert_pos]-y_object)
            
            a_coefs[0]=-al_tmp/(al_tmp+1)
            a_coefs[1]=(al_tmp-1)/al_tmp
            a_coefs[2]=1/(al_tmp**2+al_tmp)
    
            a_coefs=a_coefs/( abs(Yc[vert_pos]-y_object) )
            
            d1f[vert_pos]=(a_coefs[0]*f_wall + a_coefs[1]*f[vert_pos]+ a_coefs[2]*f[vert_pos+shift])
            
    return(d1f)


def D2_1Dy_IBM(f, settings, Yc, Y):
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
	# That means we compute the derivative of a 1D array

	# dimensions of array
    n2  = f.shape[0]
    
    d2f = np.zeros([n2])
    d1f = np.zeros([n2])
    
    d1f = np.gradient(f, Yc, edge_order=2, axis=1)
    d1f = D1_3Dy_IBM(d1f, f, settings, Yc)
    
    d2f = np.gradient( d1f, Yc, edge_order=2, axis=1)
    d2f = D1_3Dy_IBM(d2f, d1f, settings, Yc)
        
    dx2 = 2/n2
    beta = ( (Yc[0]-Y[0])**2 ) / dx2**2
    alpha = (Yc[1]-Yc[0]) / (Yc[0]-Y[0])
    a1_d = 2 / ((alpha+1)*beta)
    a2_d = -2 / (alpha*beta)
    a3_d = 2 / ((alpha**2+alpha)*beta)
    
    d2f[0] = ( a1_d*0 + a2_d*f[0] + a3_d*f[1] ) / dx2**2
        
    return(d2f)
    