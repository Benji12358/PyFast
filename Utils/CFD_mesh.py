"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Class Mesh

@description: From the settings of the current simulation data, build the mesh and get key variables (nx,ny,nz,Xc,Yc,...)
            X : Streamwise direction
            Y : Vertical direction
            Z : Spanwise direction
"""

import numpy as np
import h5py
from Settings import Settings

class CFD_mesh:
    
    def __init__(self, settings:Settings):
        """Initialisation of the CFD_mesh object.
        
        Parameters
        ----------
        settings : Settings
            User information enabling building the mesh for the simulation post-processing
    
        Returns
        -------
        CFD_mesh
            An object containing all information of the mesh.
        """

        # get a random field to initialize the mesh key points
        mean_folder = settings.root + '/' + settings.mean_path + '/' + settings.mean_path + '/Results/3D/field' \
                    + str(settings.preliminary_iterations*1000)
    
        # U (streamwise velocity here)
        fu = h5py.File(mean_folder + '/W.h5', 'r')
        U_mean = np.array(fu['W'])[:-1,:-1,:-1]
        fu.close()

        self.nx = U_mean.shape[2]
        self.ny = U_mean.shape[1]
        self.nz = U_mean.shape[0]

        del U_mean
    
        if (settings.current_path in ['PL_Vortices_LargeEps_fine', 'PL_Vortices_LargeEps_Ka']):
            self.nx = self.nx*2
            self.nz = self.nz*2

        self.X=np.linspace(0,settings.Lx,self.nx+1) # spanwise coordinates          
        self.Z=np.linspace(0,settings.Lz,self.nz+1) # streamwise coordinates
         
        self.xstep = settings.Lx/self.nx
        self.zstep = settings.Lz/self.nz

        self.Xc = np.linspace(self.xstep/2,settings.Lx-self.xstep/2,self.nx)
        self.Zc = np.linspace(self.zstep/2,settings.Lz-self.zstep/2,self.nz)

        # Wall-normal coordinates
        file = settings.root + '/' + settings.mean_path + '/' + settings.mean_path + '/Log/Mesh_generator/Y/ymesh.out'
        self.Y = np.zeros((self.ny+1))
        self.Yc = np.zeros((self.ny))

        mesh_file = open(file, 'r')
        lines = mesh_file.readlines()
        for i in range(self.ny):
            line = lines[i+1]
            array = line.split(' ')
            new_array = [x for x in array if x!='']
            self.Y[i] = float(new_array[1])
            self.Yc[i] = float(new_array[3])
        self.Y[self.ny] = 2

        ##################
        self.XY,self.YX = np.meshgrid(self.Xc,self.Yc)
        self.XZ,self.ZX = np.meshgrid(self.Xc,self.Zc-self.Zc[self.nz//2])
        self.ZY,self.YZ = np.meshgrid(self.Zc-self.Zc[self.nz//2],self.Yc)