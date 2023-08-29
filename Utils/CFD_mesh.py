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

from utilities import find_nearest, read_h5_file

class CFD_mesh:
    
    def __init__(self, settings:Settings, folder_path='3D'):
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
        mean_folder = settings.root + '/' + settings.current_path + '/' + settings.current_path + '/Results/' + str(folder_path) + '/field' \
                    + str(settings.number_iteration*1000)
                        
        # U (streamwise velocity here)
        U_mean = read_h5_file(mean_folder, 'W', settings.streamwise)

        self.nx = U_mean.shape[2]
        self.ny = U_mean.shape[1]
        self.nz = U_mean.shape[0]

        del U_mean

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
        shift = settings.shift
        settings.shift = find_nearest(self.Xc, shift)
        # self.Xc = self.Xc + shift

        self.XY,self.YX = np.meshgrid(self.Xc,self.Yc)
        self.XZ,self.ZX = np.meshgrid(self.Zc-settings.Lz/2,self.Xc)
        self.ZY,self.YZ = np.meshgrid(self.Zc-settings.Lz/2,self.Yc)
        
        ##################
        self.length_x = np.arange(1,np.floor(self.nx/2),dtype="int")
        # self.length_x = np.arange(1,2*np.floor(self.nx/4)+1,dtype="int")
        self.width_x = self.nx
        
        
    def update_mesh(self, settings:Settings):
        
        if (settings.total_domain==False):
            
            self.xs = find_nearest(self.Xc, settings.x_bounds[0])
            self.xe = find_nearest(self.Xc, settings.x_bounds[1])
            
            self.ys = find_nearest(self.Yc, settings.y_bounds[0])
            self.ye = find_nearest(self.Yc, settings.y_bounds[1])
            
            self.zs = find_nearest(self.Zc, settings.z_bounds[0])
            self.ze = find_nearest(self.Zc, settings.z_bounds[1])
            
        else:
            
            self.xs = 0
            self.xe = self.nx
            
            self.ys = 0
            self.ye = self.ny
            
            self.zs = 0
            self.ze = self.nz
            
        self.size_x = self.xe - self.xs
        self.size_y = self.ye - self.ys
        self.size_z = self.ze - self.zs