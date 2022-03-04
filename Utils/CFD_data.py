"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Class CFD_Data

@description: General class for one variable of the flow field.
            Contains mean, instantaneous and fluctuation fields.
            The spot is used to get data of a subdomain of the general channel
"""

import numpy as np

class CFD_data:
    
    def __init__(self, nx, ny, nz):
        """Initialisation of the CFD_data object.
        
        Parameters
        ----------
        nx : Integer
            Number of points in the streamwise direction
        ny : Integer
            Number of points in the vertical direction
        nz : Integer
            Number of points in the spanwise direction
    
        Returns
        -------
        CFD_data
            An object containing main information (mean, instantaneous and mean in the spot) for one flow variable
        """
        
        self.instantaneous = np.zeros((nz,ny,nx))
        self.mean = np.zeros((nz,ny,nx))
        self.mean_spot = np.zeros((ny))
        self.nz = nz
        self.ny = ny
        self.nx = nx
        
    def set_instantaneous(self, instantaneous):
        """Setting the instantaneous field of a CFD_data object.
        
        Parameters
        ----------
        instantaneous : 3D numpy array with shape (nz,ny,nx)
            Instantaneous field
    
            Update the instantaneous field of a CFD_data object
        """
        
        self.instantaneous = instantaneous
                
    def set_mean(self, mean):
        """Setting the mean field of a CFD_data object.
        
        Parameters
        ----------
        mean : 3D numpy array with shape (nz,ny,nx)
            Mean field
    
            Update the mean field of a CFD_data object
        """
        
        self.mean = mean
        
    def get_fluctuations(self):
        """Getting the fluctuated field of a CFD_data object.
    
        Returns
        -------
        3D numpy array with shape (nz,ny,nx)
            The fluctuated field of a CFD_data object from its instantaneous and mean part
        """
        
        return self.instantaneous - self.mean
    
    def set_mean_spot(self, mean_spot):
        """Setting the mean field in the spot of a CFD_data object.
        
        Parameters
        ----------
        mean : 1D numpy array with shape (ny)
            Mean field in the spot
    
            Update the mean field in the spot of a CFD_data object.
        """
        
        self.mean_spot = mean_spot
        
    def get_fluctuations_spot(self):
        """Getting the fluctuated field (compared to the main of the spot) of a CFD_data object.
    
        Returns
        -------
        3D numpy array with shape (nz,ny,nx)
            The fluctuated field of a CFD_data object from its instantaneous and mean part in the spot
        """
        
        return np.swapaxes(np.asarray([self.instantaneous[:,j,:] - self.mean_spot[j] for j in range(self.ny)]),0,1)