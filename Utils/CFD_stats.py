"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Class CFD_stats

@description: Main class of the post-processing tool.
            Use all other classes (CFD_data, CFD_mesh and Settings).
            Enable computing means, and fluctuating part of statistics.
            Enable saving statistics.
"""

import numpy as np
import h5py
from scipy.ndimage.filters import gaussian_filter
import timeit
import os

from Settings import Settings
from CFD_mesh import CFD_mesh
from CFD_data import CFD_data

from stats_utils import compute_correlations_spot, compute_correlations_x, compute_correlations_spot_3D, average_xz_spot, write_stats

from interpolation_schemes import D0s_DRP5_3Dx as interpolate_to_cc_x
from interpolation_schemes import D0s_DRP5_3Dy as interpolate_to_cc_y
from interpolation_schemes import D0s_DRP5_3Dz as interpolate_to_cc_z

from derivative_schemes import D1_3Dy, D2_3Dy
from derivative_schemes import D1_3Dx_IBM, D1_3Dy_IBM, D1_3Dz_IBM
from derivative_schemes import D1_1Dy_IBM, D2_1Dy_IBM
from derivative_schemes import D2_3Dy_IBM

from utilities import createPostProcessingDirectory, find_nearest, read_h5_file
from stats_utils import write_hdf5, write_xdmf

from plot_utils import plot_2D_contours

class CFD_stats:

    def __init__(self, settings:Settings, mesh:CFD_mesh):
        """Initialisation of the CFD_stats object.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
        Returns
        -------
        CFD_stats
            An object containing main variables (U,V,W,P and T) of a simulation
        """
        
        self.U = CFD_data(mesh.nx,mesh.ny,mesh.nz)
        self.V = CFD_data(mesh.nx,mesh.ny,mesh.nz)
        self.W = CFD_data(mesh.nx,mesh.ny,mesh.nz)
        self.P = CFD_data(mesh.nx,mesh.ny,mesh.nz)
        
        if (settings.get_T):
            self.T = CFD_data(mesh.nx,mesh.ny,mesh.nz)

    def init_ui_stats_arrays(self, settings:Settings, mesh:CFD_mesh):
        """Initialisation of the velocities stats arrays.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Update the CFD_stats object
        """
        
        # U_fluc*U_fluc in the spot
        self.UfUf_spot = np.zeros((mesh.ny))

        # V_fluc*V_fluc in the spot
        self.VfVf_spot = np.zeros((mesh.ny))

        # W_fluc*W_fluc in the spot
        self.WfWf_spot = np.zeros((mesh.ny))

        # U_fluc*V_fluc in the spot
        self.UfVf_spot = np.zeros((mesh.ny))

        # U_fluc*W_fluc in the spot
        self.UfWf_spot = np.zeros((mesh.ny))

        # V_fluc*W_fluc in the spot
        self.VfWf_spot = np.zeros((mesh.ny))

        # P_fluc*U_fluc in the spot
        self.PfUf_spot = np.zeros((mesh.ny))

        ####################################################
        # U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUf_spot = np.zeros((mesh.ny))

        # U_fluc*U_fluc*V_fluc in the spot
        self.UfUfVf_spot = np.zeros((mesh.ny))

        # U_fluc*V_fluc*V_fluc in the spot
        self.UfVfVf_spot = np.zeros((mesh.ny))

        # V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVf_spot = np.zeros((mesh.ny))

        # W_fluc*W_fluc*V_fluc in the spot
        self.WfWfVf_spot = np.zeros((mesh.ny))

        # W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWf_spot = np.zeros((mesh.ny))

        ####################################################
        # U_fluc*U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUfUf_spot = np.zeros((mesh.ny))

        # V_fluc*V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVfVf_spot = np.zeros((mesh.ny))

        # W_fluc*W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWfWf_spot = np.zeros((mesh.ny))
            
        if (settings.compute_vorticity):
            
            self.Vort_fluc = {
                'name': 'Vort_fluc',
                'Wall normal coordinate': mesh.Yc,
                'omega_x': np.zeros((mesh.ny)),
                'omega_y': np.zeros((mesh.ny)),
                'omega_z': np.zeros((mesh.ny)),
                'omega_x rms': np.zeros((mesh.ny)),
                'omega_y rms': np.zeros((mesh.ny)),
                'omega_z rms': np.zeros((mesh.ny))
                }
        
        if (settings.compute_transport_equation_terms):
            
            self.uu_budget = {
            'name': 'uu_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
            self.vv_budget = {
            'name': 'vv_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
            self.ww_budget = {
            'name': 'ww_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
            self.uv_budget = {
            'name': 'uv_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
        if (settings.compute_spectrums):
            
            self.uu_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))
            self.vv_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))
            self.ww_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))
            self.uv_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))

    def init_local_stats_arrays(self, settings:Settings, mesh:CFD_mesh):
        """Initialisation of the local stats arrays.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Update the CFD_stats object
        """
        
        # # Re_tau(x)
        # self.Re_tau_x = np.zeros((mesh.nx))

        # # T_tau(x)
        # self.T_tau_x = np.zeros((mesh.nx))

        # # Nu(x)
        # self.Nu_x = np.zeros((mesh.nx))

        # Ek(x)
        self.Ek_x = np.zeros((mesh.nx))

        # TT(x)
        self.TT_x = np.zeros((mesh.nx))

        # Omega(x)
        self.Omega_x = np.zeros((mesh.nx))

    def init_ui_stats_arrays_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Initialisation of the velocities stats arrays.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Update the CFD_stats object
        """
        
        # U_fluc*U_fluc in the spot
        self.UfUf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # V_fluc*V_fluc in the spot
        self.VfVf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # W_fluc*W_fluc in the spot
        self.WfWf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # U_fluc*V_fluc in the spot
        self.UfVf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # U_fluc*W_fluc in the spot
        self.UfWf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # V_fluc*W_fluc in the spot
        self.VfWf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # P_fluc*U_fluc in the spot
        self.PfUf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        ####################################################
        # U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # U_fluc*U_fluc*V_fluc in the spot
        self.UfUfVf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # U_fluc*V_fluc*V_fluc in the spot
        self.UfVfVf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # W_fluc*W_fluc*V_fluc in the spot
        self.WfWfVf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        ####################################################
        # U_fluc*U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUfUf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # V_fluc*V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVfVf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))

        # W_fluc*W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWfWf_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
        if (settings.compute_vorticity):
            
            self.Vort_fluc = {
                'name': 'Vort_fluc',
                'Wall normal coordinate': mesh.Yc,
                'omega_x': np.zeros((mesh.ny)),
                'omega_y': np.zeros((mesh.ny)),
                'omega_z': np.zeros((mesh.ny)),
                'omega_x rms': np.zeros((mesh.ny)),
                'omega_y rms': np.zeros((mesh.ny)),
                'omega_z rms': np.zeros((mesh.ny))
                }
        
        if (settings.compute_transport_equation_terms):
            
            self.uu_budget = {
            'name': 'uu_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
            self.vv_budget = {
            'name': 'vv_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
            self.ww_budget = {
            'name': 'ww_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
            self.uv_budget = {
            'name': 'uv_budget',
            'Wall normal coordinate': mesh.Yc,
            'Production': np.zeros((mesh.ny)),
            'Turbulent Diffusive Flux': np.zeros((mesh.ny)),
            'Pressure velocity gradient correlation': np.zeros((mesh.ny)),
            'Molecular diffusion': np.zeros((mesh.ny)),
            'Dissipation': np.zeros((mesh.ny))
            }
            
        if (settings.compute_spectrums):
            
            self.uu_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))
            self.vv_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))
            self.ww_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))
            self.uv_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))
            
    
    def init_T_stats_arrays(self, settings:Settings, mesh:CFD_mesh):
        """Initialisation of the temperature stats arrays.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Update the CFD_stats object
        """
        
        ####################################################
        # T_fluc*T_fluc in the spot
        self.TfTf_spot = np.zeros((mesh.ny))

        # U_fluc*T_fluc in the spot
        self.UfTf_spot = np.zeros((mesh.ny))

        # V_fluc*T_fluc in the spot
        self.VfTf_spot = np.zeros((mesh.ny))

        # W_fluc*T_fluc in the spot
        self.WfTf_spot = np.zeros((mesh.ny))
        
        ####################################################
        # T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTf_spot = np.zeros((mesh.ny))
        
        # V_fluc*T_fluc*T_fluc in the spot
        self.VfTfTf_spot = np.zeros((mesh.ny))

        ####################################################
        # T_fluc*T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTfTf_spot = np.zeros((mesh.ny))
        
        if (settings.compute_temperature_budgets):
            
            self.tt_budget = {
                'name': 'TemperatureBudgetsSpot',
                'Wall normal coordinate': mesh.Yc,
                'Advection': np.zeros((mesh.ny)),
                'Turbulent Transport': np.zeros((mesh.ny)),
                'Production': np.zeros((mesh.ny)),
                'Molecular diffusion': np.zeros((mesh.ny)),
                'Dissipation': np.zeros((mesh.ny))
                }
                
        if (settings.compute_spectrums):
            
            self.tt_spectrum = np.zeros((mesh.nz//2-1,mesh.ny//2,mesh.nx//2-1))


    def average_ui_stats_arrays(self, settings:Settings):
        """Average the velocities stats arrays by the total number of iterations.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
    
            Update the CFD_stats object
        """
        
        # U_fluc*U_fluc in the spot
        self.UfUf_spot = self.UfUf_spot/len(settings.all_iterations)

        # V_fluc*V_fluc in the spot
        self.VfVf_spot = self.VfVf_spot/len(settings.all_iterations)

        # W_fluc*W_fluc in the spot
        self.WfWf_spot = self.WfWf_spot/len(settings.all_iterations)

        # U_fluc*V_fluc in the spot
        self.UfVf_spot = self.UfVf_spot/len(settings.all_iterations)

        # U_fluc*W_fluc in the spot
        self.UfWf_spot = self.UfWf_spot/len(settings.all_iterations)

        # V_fluc*W_fluc in the spot
        self.VfWf_spot = self.VfWf_spot/len(settings.all_iterations)

        # P_fluc*U_fluc in the spot
        self.PfUf_spot = self.PfUf_spot/len(settings.all_iterations)

        ####################################################
        # U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUf_spot = self.UfUfUf_spot/len(settings.all_iterations)

        # U_fluc*U_fluc*V_fluc in the spot
        self.UfUfVf_spot = self.UfUfVf_spot/len(settings.all_iterations)

        # U_fluc*V_fluc*V_fluc in the spot
        self.UfVfVf_spot = self.UfVfVf_spot/len(settings.all_iterations)

        # V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVf_spot = self.VfVfVf_spot/len(settings.all_iterations)

        # W_fluc*W_fluc*V_fluc in the spot
        self.WfWfVf_spot = self.WfWfVf_spot/len(settings.all_iterations)

        # W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWf_spot = self.WfWfWf_spot/len(settings.all_iterations)

        ####################################################
        # U_fluc*U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUfUf_spot = self.UfUfUfUf_spot/len(settings.all_iterations)

        # V_fluc*V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVfVf_spot = self.VfVfVfVf_spot/len(settings.all_iterations)

        # W_fluc*W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWfWf_spot = self.WfWfWfWf_spot/len(settings.all_iterations)
            
        if (settings.compute_vorticity):
            
            self.Vort_fluc['omega_x'] = self.Vort_fluc['omega_x']/len(settings.all_iterations)
            self.Vort_fluc['omega_y'] = self.Vort_fluc['omega_y']/len(settings.all_iterations)
            self.Vort_fluc['omega_z'] = self.Vort_fluc['omega_z']/len(settings.all_iterations)
            self.Vort_fluc['omega_x rms'] = self.Vort_fluc['omega_x rms']/len(settings.all_iterations)
            self.Vort_fluc['omega_y rms'] = self.Vort_fluc['omega_y rms']/len(settings.all_iterations)
            self.Vort_fluc['omega_z rms'] = self.Vort_fluc['omega_z rms']/len(settings.all_iterations)
        
        if (settings.compute_transport_equation_terms):
            
            self.uu_budget['Production'] = self.uu_budget['Production']/len(settings.all_iterations)
            self.uu_budget['Turbulent Diffusive Flux'] = self.uu_budget['Turbulent Diffusive Flux']/len(settings.all_iterations)
            self.uu_budget['Pressure velocity gradient correlation'] = self.uu_budget['Pressure velocity gradient correlation']/len(settings.all_iterations)
            self.uu_budget['Molecular diffusion'] = self.uu_budget['Molecular diffusion']/len(settings.all_iterations)
            self.uu_budget['Dissipation'] = self.uu_budget['Dissipation']/len(settings.all_iterations)
            
            self.vv_budget['Production'] = self.vv_budget['Production']/len(settings.all_iterations)
            self.vv_budget['Turbulent Diffusive Flux'] = self.vv_budget['Turbulent Diffusive Flux']/len(settings.all_iterations)
            self.vv_budget['Pressure velocity gradient correlation'] = self.vv_budget['Pressure velocity gradient correlation']/len(settings.all_iterations)
            self.vv_budget['Molecular diffusion'] = self.vv_budget['Molecular diffusion']/len(settings.all_iterations)
            self.vv_budget['Dissipation'] = self.vv_budget['Dissipation']/len(settings.all_iterations)
            
            self.ww_budget['Production'] = self.ww_budget['Production']/len(settings.all_iterations)
            self.ww_budget['Turbulent Diffusive Flux'] = self.ww_budget['Turbulent Diffusive Flux']/len(settings.all_iterations)
            self.ww_budget['Pressure velocity gradient correlation'] = self.ww_budget['Pressure velocity gradient correlation']/len(settings.all_iterations)
            self.ww_budget['Molecular diffusion'] = self.ww_budget['Molecular diffusion']/len(settings.all_iterations)
            self.ww_budget['Dissipation'] = self.ww_budget['Dissipation']/len(settings.all_iterations)
            
            self.uv_budget['Production'] = self.uv_budget['Production']/len(settings.all_iterations)
            self.uv_budget['Turbulent Diffusive Flux'] = self.uv_budget['Turbulent Diffusive Flux']/len(settings.all_iterations)
            self.uv_budget['Pressure velocity gradient correlation'] = self.uv_budget['Pressure velocity gradient correlation']/len(settings.all_iterations)
            self.uv_budget['Molecular diffusion'] = self.uv_budget['Molecular diffusion']/len(settings.all_iterations)
            self.uv_budget['Dissipation'] = self.uv_budget['Dissipation']/len(settings.all_iterations)
            
        if (settings.compute_spectrums):
            
            self.uu_spectrum = self.uu_spectrum/len(settings.all_iterations)
            self.vv_spectrum = self.vv_spectrum/len(settings.all_iterations)
            self.ww_spectrum = self.ww_spectrum/len(settings.all_iterations)
            self.uv_spectrum = self.uv_spectrum/len(settings.all_iterations)


    def average_local_stats_arrays(self, settings:Settings):
        """Average the local stats arrays by the total number of iterations.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
    
            Update the CFD_stats object
        """
        
        # # Re_tau(x)
        # self.Re_tau_x = self.Re_tau_x/len(settings.all_iterations)

        # # T_tau(x)
        # self.T_tau_x = self.T_tau_x/len(settings.all_iterations)

        # # Nu(x)
        # self.Nu_x = self.Nu_x/len(settings.all_iterations)

        # Ek(x)
        self.Ek_x = self.Ek_x/len(settings.all_iterations)

        # TT(x)
        self.TT_x = self.TT_x/len(settings.all_iterations)

        # Omega(x)
        self.Omega_x = self.Omega_x/len(settings.all_iterations)


    def average_T_stats_arrays(self, settings:Settings):
        """Average the temperature stats arrays by the total number of iterations.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
    
            Update the CFD_stats object
        """
        
        ####################################################
        # T_fluc*T_fluc in the spot
        self.TfTf_spot = self.TfTf_spot/len(settings.all_iterations)

        # U_fluc*T_fluc in the spot
        self.UfTf_spot = self.UfTf_spot/len(settings.all_iterations)

        # V_fluc*T_fluc in the spot
        self.VfTf_spot = self.VfTf_spot/len(settings.all_iterations)

        # W_fluc*T_fluc in the spot
        self.WfTf_spot = self.WfTf_spot/len(settings.all_iterations)
        
        ####################################################
        # T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTf_spot = self.TfTfTf_spot/len(settings.all_iterations)
        
        # V_fluc*T_fluc*T_fluc in the spot
        self.VfTfTf_spot = self.VfTfTf_spot/len(settings.all_iterations)

        ####################################################
        # T_fluc*T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTfTf_spot = self.TfTfTfTf_spot/len(settings.all_iterations)
        
        if (settings.compute_temperature_budgets):
            
            self.tt_budget['Advection'] = self.tt_budget['Advection']/len(settings.all_iterations)
            self.tt_budget['Turbulent Transport'] = self.tt_budget['Turbulent Transport']/len(settings.all_iterations)
            self.tt_budget['Production'] = self.tt_budget['Production']/len(settings.all_iterations)
            self.tt_budget['Molecular diffusion'] = self.tt_budget['Molecular diffusion']/len(settings.all_iterations)
            self.tt_budget['Dissipation'] = self.tt_budget['Dissipation']/len(settings.all_iterations)
            
        if (settings.compute_spectrums):
            
            self.tt_spectrum = self.tt_spectrum/len(settings.all_iterations)


    def shift_fields_in_x(self, settings:Settings, mesh:CFD_mesh):
        """Shift instantaneous fields in x direction
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
        Returns
        -------
        CFD_stats
            An object containing main variables (U,V,W,P and T) of a simulation
        """
        
        self.U.set_instantaneous(np.roll(self.U.instantaneous,-settings.shift,axis=-1))
        self.V.set_instantaneous(np.roll(self.V.instantaneous,-settings.shift,axis=-1))
        self.W.set_instantaneous(np.roll(self.W.instantaneous,-settings.shift,axis=-1))
        self.P.set_instantaneous(np.roll(self.P.instantaneous,-settings.shift,axis=-1))

        self.spot = np.roll(self.spot,-settings.shift,axis=-1)
        
        if (settings.get_T):
            self.T.set_instantaneous(np.roll(self.T.instantaneous,-settings.shift,axis=-1))

            self.thermal_spot = np.roll(self.thermal_spot,-settings.shift,axis=-1)
        
        

    def set_mean_fields(self, settings:Settings, mesh:CFD_mesh):
        """Given the path of the mean simulation, set the mean fields of all variables.
        Only used for transitional flows.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Set the mean fields of the CFD_stats object
        """
        
        mean_folder = settings.root + '/' + settings.mean_path + '/' + settings.mean_path + '/Results/3D/field' + str(settings.preliminary_iterations*1000)

        # U (streamwise velocity here)
        U_mean = read_h5_file(mean_folder, 'W', settings.streamwise)
        # V
        V_mean = read_h5_file(mean_folder, 'V', settings.streamwise)
        # U (spanwise velocity here)
        W_mean = read_h5_file(mean_folder, 'U', settings.streamwise)
        # Pressure
        P_mean = read_h5_file(mean_folder, 'P', settings.streamwise)
        
        # from staggered to cell center field
        U_mean = interpolate_to_cc_x(U_mean)
        V_mean = interpolate_to_cc_y(V_mean)
        W_mean = interpolate_to_cc_z(W_mean)
        
        if (settings.get_T):
            T_mean = read_h5_file(mean_folder, 'sca1', settings.streamwise)
        
        print(mean_folder)
        
        if (settings.current_path in ['PL_Vortices_LargeEps_fine', 'New_PL_Vort', 'PL_Vortices_LargeEps_Ka', 'PL_Vort_LargeEps_big']):
            
            U_backup = U_mean[mesh.nz//8-1,:,mesh.nx//8-1]
            P_backup = P_mean[mesh.nz//8-1,:,mesh.nx//8-1]
                                            
            W_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            V_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            U_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            P_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
            for i in range(mesh.nx):
                for k in range(mesh.nz):
                    U_mean[k,:,i] = U_backup
                    P_mean[k,:,i] = P_backup
                    
            if (settings.get_T):
                
                T_backup = T_mean[mesh.nz//8-1,:,mesh.nx//8-1]
                
                T_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
                for i in range(mesh.nx):
                    for k in range(mesh.nz):
                        T_mean[k,:,i] = T_backup
                        
        self.U.set_mean(U_mean[:,:,:mesh.nx])
        self.V.set_mean(V_mean[:,:,:mesh.nx])
        self.W.set_mean(W_mean[:,:,:mesh.nx])
        self.P.set_mean(P_mean[:,:,:mesh.nx])
            
        del U_mean, V_mean, W_mean, P_mean
            
        if ('Ka' in settings.current_path) and (settings.get_T):
                
            T_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        
            for j in range(mesh.ny):
                    T_mean[:,j,:] = settings.Pr * np.sqrt(2*settings.Re) * (- 0.5*mesh.Yc[j]**3 + 1/8*mesh.Yc[j]**4 + mesh.Yc[j])
        
        if (settings.get_T):
            
            self.T.set_mean(T_mean[:,:,:mesh.nx])
                
            del T_mean
            
            
    def compute_mean_fields(self, settings:Settings, mesh:CFD_mesh, folder_path='3D'):
        """Compute the mean fields of all variables by using a time-averaging procedure.
        Can be used for all kind of flows.
        Compute also the mean in the spot (1D numpy array with shape (ny)) and basic statistics numbers (Re_tau, T_tau, Nu)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Compute the mean fields of the CFD_stats object
        """
        
        print(settings.root + '/' + settings.current_path + '/' + settings.current_path)
        print('Averaging from ' + str(settings.all_iterations[0]) + ' to ' + str(settings.all_iterations[-1]))
        
        U_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        V_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        W_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        P_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        
        if (settings.get_T):
            T_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        
        for i in settings.all_iterations:
            
            mean_folder = settings.root + '/' + settings.current_path + '/' + settings.current_path + '/Results/' + str(folder_path) + '/field' + str(i*1000)
    
            # U (streamwise velocity here)
            U_tmp = read_h5_file(mean_folder, 'W', settings.streamwise)
            # V
            V_tmp = read_h5_file(mean_folder, 'V', settings.streamwise)
            # U (spanwise velocity here)
            W_tmp = read_h5_file(mean_folder, 'U', settings.streamwise)
            # Pressure
            P_mean += read_h5_file(mean_folder, 'P', settings.streamwise)
            
            # from staggered to cell center field
            if (settings.IBM_flag):
                U_mean += interpolate_to_cc_x(U_tmp*self.spot)[:,:,:mesh.nx]
                V_mean += interpolate_to_cc_y(V_tmp*self.spot)[:,:,:mesh.nx]
                W_mean += interpolate_to_cc_z(W_tmp*self.spot)[:,:,:mesh.nx]
            else:
                U_mean += interpolate_to_cc_x(U_tmp)[:,:,:mesh.nx]
                V_mean += interpolate_to_cc_y(V_tmp)[:,:,:mesh.nx]
                W_mean += interpolate_to_cc_z(W_tmp)[:,:,:mesh.nx]
            
            if (settings.get_T):
                T_mean += read_h5_file(mean_folder, 'sca1', settings.streamwise)
        
        self.U.set_mean(U_mean/len(settings.all_iterations))
        self.V.set_mean(V_mean/len(settings.all_iterations))
        self.W.set_mean(W_mean/len(settings.all_iterations))
        self.P.set_mean(P_mean[:,:,:mesh.nx]/len(settings.all_iterations))
            
        del U_mean, V_mean, W_mean, P_mean

        ####################################################
        # U_mean in the spot
        tmp = self.spot * self.U.mean
        tmp[tmp==0] = np.nan
        U_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        # V_mean in the spot
        tmp = self.spot * self.V.mean
        tmp[tmp==0] = np.nan
        V_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        # W_mean in the spot
        tmp = self.spot * self.W.mean
        tmp[tmp==0] = np.nan
        W_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        # P_mean in the spot
        tmp = self.spot * self.P.mean
        tmp[tmp==0] = np.nan
        P_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        self.U.set_mean_spot(U_mean_spot)
        self.V.set_mean_spot(V_mean_spot)
        self.W.set_mean_spot(W_mean_spot)
        self.P.set_mean_spot(P_mean_spot)
        
        # Re_tau in the spot
        self.Re_tau_spot = np.sqrt( settings.Re * self.U.mean_spot[0] / mesh.Yc[0])
        
        if (settings.symmetry):
            self.Re_tau_spot += np.sqrt( settings.Re * self.U.mean_spot[-1] / mesh.Yc[0])
            self.Re_tau_spot = self.Re_tau_spot/2
        
        del U_mean_spot, V_mean_spot, W_mean_spot, P_mean_spot
        
        if (settings.get_T):
                            
            self.T.set_mean(T_mean/len(settings.all_iterations))
                
            del T_mean
            
            ####################################################
            # T_mean in the spot
            tmp = self.thermal_spot * self.T.mean
            tmp[tmp==0] = np.nan
            T_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
    
            self.T.set_mean_spot(T_mean_spot)
            
            del T_mean_spot
            
            ####################################################
            # T_tau in the spot
            self.T_tau_spot = abs((self.T.mean_spot[0]+settings.deltaT) / (mesh.Yc[0]*self.Re_tau_spot*settings.Pr))
        
            if (settings.symmetry):
                self.T_tau_spot += abs((self.T.mean_spot[-1]-settings.deltaT) / (mesh.Yc[0]*self.Re_tau_spot*settings.Pr))
                self.T_tau_spot = self.T_tau_spot/2
                    
            #### Now for Nusselt
            # u_theta = 0
            # u = 0
            
            # for k in range(mesh.nz):
            #     for j in range(mesh.ny//2):
                    
            #         u_theta += self.thermal_spot[k,j,-1]*self.U.mean[k,j,-1]*(self.T.mean[k,j,-1]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            #         u += self.thermal_spot[k,j,-1]*self.U.mean[k,j,-1]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            # theta_b = (u_theta/u) * (1/self.T_tau_spot)
            
            # self.Nu_spot = 2*self.Re_tau_spot*settings.Pr/theta_b
            
            u_theta = 0
            u = 0
            
            for j in range(mesh.ny//2):
                
                u_theta += self.U.mean_spot[j]*(self.T.mean_spot[j]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])
                
                u += self.U.mean_spot[j]*(mesh.Y[j+1]-mesh.Y[j])
                    
            theta_b = (u_theta/u) * (1/self.T_tau_spot)
            
            self.Nu_spot = 2*self.Re_tau_spot*settings.Pr/theta_b
                
            if ('PSK' in settings.current_path) or ('Ka' in settings.current_path):
                # Kawamura temperature
                self.T_tau_spot = 1 
                    
                #### Now for Nusselt
                u_theta = 0
                u = 0
                
                # for k in range(mesh.nz):
                #     for j in range(mesh.ny):
                        
                #         u_theta += self.thermal_spot[k,j,-1]*self.U.mean[k,j,-1]*self.T.mean[k,j,-1]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                        
                #         u += self.thermal_spot[k,j,-1]*self.U.mean[k,j,-1]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                
                for k in range(mesh.nz):
                    for j in range(mesh.ny//2):
                        
                        u_theta += self.U.mean_spot[j]*self.T.mean_spot[j]*(mesh.Y[j+1]-mesh.Y[j])
                        
                        u += self.U.mean_spot[j]*(mesh.Y[j+1]-mesh.Y[j])
                        
                # print(u_theta, u)
                theta_b = (u_theta/u) * (1/self.T_tau_spot)
                
                self.Nu_spot = 2*self.Re_tau_spot*settings.Pr/theta_b
                
            # # Nusselt in the spot
            # tmp = self.thermal_spot[:,0,:] * self.T.mean[:,0,:]
            # tmp[tmp==0] = np.nan
            # self.Nu_spot = abs((np.nanmean(tmp)+settings.deltaT) / (mesh.Yc[0])) * ( 1 / settings.deltaT )
        
            # if (settings.symmetry):
            #     tmp = self.thermal_spot[:,-1,:] * self.T.mean[:,-1,:]
            #     tmp[tmp==0] = np.nan
            #     self.Nu_spot += abs((np.nanmean(tmp)-settings.deltaT) / (mesh.Yc[0])) * ( 1 / settings.deltaT )
            #     self.Nu_spot = self.Nu_spot/2
                    
                
            
    def compute_drag_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute the mean fields of all variables by using a time-averaging procedure.
        Can be used for all kind of flows.
        Compute also the mean in the spot (1D numpy array with shape (ny)) and basic statistics numbers (Re_tau, T_tau, Nu)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Compute the mean fields of the CFD_stats object
        """
        
        friction_lw = 0
        friction_rc = 0

        for i in range(mesh.nx):
            for k in range(mesh.nz):

                friction_lw += ( self.spot[k,0,i] * self.U.mean[k,0,i]/mesh.Yc[0] ) * mesh.xstep * mesh.zstep

        # for both roughness on lower wall TOP
        for n in range(len(settings.i_start)):

            j_s  = settings.j_start[n]
            j_e  = settings.j_end[n]
            y_end = settings.y_end[n]

            i_s  = settings.i_start[n]
            i_e  = settings.i_end[n]

            k_s  = settings.k_start[n]
            k_e  = settings.k_end[n]

            if (j_s==0):

                tmp_friction = 0

                for i in range(i_s,min(i_e+1,mesh.nx)):
                    for k in range(k_s,min(k_e+1,mesh.nz)):

                        tmp_friction += ( (1-self.spot[k,j_e,i]) * self.U.mean[k,j_e+1,i]/(mesh.Yc[j_e+1]-y_end)) * mesh.xstep * mesh.zstep

                friction_rc += tmp_friction
        
        surface = settings.Lz*settings.Lx
        volume_fluide = 2 * surface - len(settings.i_start) * settings.y_end[0]**2 * settings.z_end[0]
        
        self.friction_drag = 2 * ( abs(friction_lw) + friction_rc ) / ( settings.Re * surface )
        self.total_drag = ( (self.Re_tau_spot/settings.Re)**2 ) * ( volume_fluide / surface )
        self.pressure_drag = self.total_drag - self.friction_drag
                    
                
            
    def compute_T_tau_and_Nu_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute the mean fields of all variables by using a time-averaging procedure.
        Can be used for all kind of flows.
        Compute also the mean in the spot (1D numpy array with shape (ny)) and basic statistics numbers (Re_tau, T_tau, Nu)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Compute the mean fields of the CFD_stats object
        """
        
        if ('PSK' in settings.current_path):
        # Kawamura temperature
            self.T_tau_spot = 1
            
            #### Now for Nusselt
            u_theta = 0
            u = 0
            
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta += self.thermal_spot[k,j,0]*self.U.mean[k,j,0]*(self.T.mean[k,j,0]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u += self.thermal_spot[k,j,0]*self.U.mean[k,j,0]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            theta_b = (u_theta/u)
            
            self.Nu_spot = 2*self.Re_tau_spot*settings.Pr/theta_b
                
        else:
            
            ########################
            ###### calcul global
            ########################
            
            print('calcul global')
                        
            heat_flux_lw = 0
            heat_flux_rc = 0
            rc_surface = 0
            heat_flux = 0
                        
            for i in range(mesh.nx):
                for k in range(mesh.nz):
                    
                    heat_flux_lw += ( self.thermal_spot[k,0,i] * (self.T.mean[k,0,i]+settings.deltaT)/mesh.Yc[0] ) * mesh.xstep * mesh.zstep
                                
            # for both roughness on lower wall TOP
            for n in range(len(settings.i_start)):
                
                j_s  = settings.j_start[n]
                j_e  = settings.j_end[n]
                y_end = settings.y_end[n]
                
                i_s  = settings.i_start[n]
                i_e  = settings.i_end[n]
                
                k_s  = settings.k_start[n]
                k_e  = settings.k_end[n]
                
                if (j_s==0):
            
                    tmp_heat_flux = 0
                    rc_surface += (settings.x_end[n]-settings.x_start[n]) * (settings.z_end[n]-settings.z_start[n])
                    
                    for i in range(i_s,min(i_e+1,mesh.nx)):
                        for k in range(k_s,min(k_e+1,mesh.nz)):
                            
                            tmp_heat_flux += ( (1-self.thermal_spot[k,j_e,i]) * (self.T.mean[k,j_e+1,i]+settings.deltaT)/(mesh.Yc[j_e+1]-y_end)) * mesh.xstep * mesh.zstep
                            
                    heat_flux_rc += tmp_heat_flux
                        
            print('rc_surface (%tot) =', rc_surface/(settings.Lz*settings.Lx))
            # print(heat_flux_lw, heat_flux_rc)
            heat_flux = ( heat_flux_lw + heat_flux_rc ) / (settings.Lz*settings.Lx)
            
            self.T_tau_spot = heat_flux/(self.Re_tau_spot*settings.Pr)
            
            print('T_tau', self.T_tau_spot)
            
            #### Now for Nusselt
            u_theta = 0
            u = 0
                    
            ## Be carefull, if Lyons, then we use only half of the channel
            for j in range(mesh.ny//2):
                    
                u_theta += self.U.mean_spot[j]*(self.T.mean_spot[j]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])
                u += self.U.mean_spot[j]*(mesh.Y[j+1]-mesh.Y[j])
                    
            self.theta_b = (u_theta/u) * (1/self.T_tau_spot)
            self.Nu_spot = 2*self.Re_tau_spot*settings.Pr/self.theta_b
            
            # ########### VERSION SEDAT
            # self.Nu_spot = heat_flux        
            
            # print('Nu', self.Nu_spot)
                    
                
            
    def compute_T_tau_and_Nu_ibm_local(self, settings:Settings, mesh:CFD_mesh, place):
        """Compute the mean fields of all variables by using a time-averaging procedure.
        Can be used for all kind of flows.
        Compute also the mean in the spot (1D numpy array with shape (ny)) and basic statistics numbers (Re_tau, T_tau, Nu)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Compute the mean fields of the CFD_stats object
        """
        
        print(place)
        if (place=='between'):
            # T_tau is already well computed
            
            # Computation of Re_tau
            # self.Re_tau_spot = np.sqrt( settings.Re * abs(self.U.mean_spot[0]) / (mesh.Yc[0]))
            
            # Computation of Nu
            u_theta = 0
            u = 0
            
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta += self.thermal_spot[k,j,-1]*self.U.mean[k,j,-1]*(self.T.mean[k,j,-1]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u += self.thermal_spot[k,j,-1]*self.U.mean[k,j,-1]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            theta_b = (u_theta/u) * (1/self.T_tau_spot)
            
            self.Nu_spot = 2*self.Re_tau_spot*settings.Pr/theta_b
            
            
        elif (place=='above'):
            
            # In this case, we need to compute correctly both Re_tau, T_tau and Nu above roughness
            j_e  = settings.j_end[0]
            y_end = settings.y_end[0]
            
            # Computation of Re_tau
            self.Re_tau_spot = np.sqrt( settings.Re * self.U.mean_spot[j_e+1] / (mesh.Yc[j_e+1]-y_end))
       
            # Computation of T_tau
            self.T_tau_spot = (self.T.mean_spot[j_e+1]+settings.deltaT)/((mesh.Yc[j_e+1]-y_end)*self.Re_tau_spot*settings.Pr)
            
            # Computation of Nu
            # self.Nu_spot = (self.T.mean_spot[j_e+1]+settings.deltaT)/(mesh.Yc[j_e+1]-y_end)
            u_theta = 0
            u = 0
            
            self.U.mean_spot[:j_e+1] = 0
                    
            ## Be carefull, if Lyons, then we use only half of the channel
            for j in range(mesh.ny//2):
                    
                u_theta += self.U.mean_spot[j]*(self.T.mean_spot[j]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])
                u += self.U.mean_spot[j]*(mesh.Y[j+1]-mesh.Y[j])
                    
            theta_b = (u_theta/u) * (1/self.T_tau_spot)
            self.Nu_spot = 2*self.Re_tau_spot*settings.Pr/theta_b
               
        else:
            
            print("Error place")
                
            
    def compute_mean_fields_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute the mean fields of all variables by using a time-averaging procedure.
        Can be used for all kind of flows.
        Compute also the mean in the spot (1D numpy array with shape (ny)) and basic statistics numbers (Re_tau, T_tau, Nu)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Compute the mean fields of the CFD_stats object
        """
        
        tmp_mask_ibm = np.ones((mesh.nz,mesh.ny,mesh.nx))
        
        for n in range(len(settings.i_start)):
            
            if settings.j_start[n]==0:
            
                tmp_mask_ibm[settings.k_start[n]:settings.k_end[n]+1,:mesh.ny//2+1,settings.i_start[n]:settings.i_end[n]+1] = 0
                
            else:
            
                tmp_mask_ibm[settings.k_start[n]:settings.k_end[n]+1,mesh.ny//2+1:mesh.ny+1,settings.i_start[n]:settings.i_end[n]+1] = 0
             
        # Out roughness      
        tmp = tmp_mask_ibm * self.U.mean
        tmp[tmp==0] = np.nan
        U_mean_wo_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        tmp = tmp_mask_ibm * self.V.mean
        tmp[tmp==0] = np.nan
        V_mean_wo_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        tmp = tmp_mask_ibm * self.W.mean
        tmp[tmp==0] = np.nan
        W_mean_wo_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        tmp = tmp_mask_ibm * self.P.mean
        tmp[tmp==0] = np.nan
        P_mean_wo_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # In roughness
        tmp = (1-tmp_mask_ibm) * self.U.mean
        tmp[tmp==0] = np.nan
        U_mean_w_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        tmp = (1-tmp_mask_ibm) * self.V.mean
        tmp[tmp==0] = np.nan
        V_mean_w_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        tmp = (1-tmp_mask_ibm) * self.W.mean
        tmp[tmp==0] = np.nan
        W_mean_w_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        tmp = (1-tmp_mask_ibm) * self.P.mean
        tmp[tmp==0] = np.nan
        P_mean_w_roughness = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Out roughness
        u3D_mean_wo_roughness = np.reshape(U_mean_wo_roughness,(mesh.ny,1))
        u3D_mean_wo_roughness = np.tile(u3D_mean_wo_roughness,(mesh.nz,1,mesh.nx))
        
        v3D_mean_wo_roughness = np.reshape(V_mean_wo_roughness,(mesh.ny,1))
        v3D_mean_wo_roughness = np.tile(v3D_mean_wo_roughness,(mesh.nz,1,mesh.nx))
        
        w3D_mean_wo_roughness = np.reshape(W_mean_wo_roughness,(mesh.ny,1))
        w3D_mean_wo_roughness = np.tile(w3D_mean_wo_roughness,(mesh.nz,1,mesh.nx))
        
        p3D_mean_wo_roughness = np.reshape(P_mean_wo_roughness,(mesh.ny,1))
        p3D_mean_wo_roughness = np.tile(p3D_mean_wo_roughness,(mesh.nz,1,mesh.nx))
        
        # In roughness
        u3D_mean_w_roughness = np.reshape(U_mean_w_roughness,(mesh.ny,1))
        u3D_mean_w_roughness = np.tile(u3D_mean_w_roughness,(mesh.nz,1,mesh.nx))
        
        v3D_mean_w_roughness = np.reshape(V_mean_w_roughness,(mesh.ny,1))
        v3D_mean_w_roughness = np.tile(v3D_mean_w_roughness,(mesh.nz,1,mesh.nx))
        
        w3D_mean_w_roughness = np.reshape(W_mean_w_roughness,(mesh.ny,1))
        w3D_mean_w_roughness = np.tile(w3D_mean_w_roughness,(mesh.nz,1,mesh.nx))
        
        p3D_mean_w_roughness = np.reshape(P_mean_w_roughness,(mesh.ny,1))
        p3D_mean_w_roughness = np.tile(p3D_mean_w_roughness,(mesh.nz,1,mesh.nx))
        
        self.U.set_mean( tmp_mask_ibm*u3D_mean_wo_roughness + (1-tmp_mask_ibm)*u3D_mean_w_roughness )
        self.V.set_mean( tmp_mask_ibm*v3D_mean_wo_roughness + (1-tmp_mask_ibm)*v3D_mean_w_roughness )
        self.W.set_mean( tmp_mask_ibm*w3D_mean_wo_roughness + (1-tmp_mask_ibm)*w3D_mean_w_roughness )
        self.P.set_mean( tmp_mask_ibm*p3D_mean_wo_roughness + (1-tmp_mask_ibm)*p3D_mean_w_roughness )
        
            
            
    def average_in_patterns(self, settings:Settings, mesh:CFD_mesh):
        """Compute the mean fields of all variables by using a time-averaging procedure.
        Can be used for all kind of flows.
        Compute also the mean in the spot (1D numpy array with shape (ny)) and basic statistics numbers (Re_tau, T_tau, Nu)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Compute the mean fields of the CFD_stats object
        """
            
        U_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        V_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        W_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        P_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        
        if (settings.get_T):
            T_mean = np.zeros((mesh.nz, mesh.ny, mesh.nx))
    
        num_of_pattern = len(settings.x_start)//4
    
        print('Starting procedure with ' + str(num_of_pattern) + ' patterns')
        
        nx_pattern = mesh.nx//num_of_pattern
            
        U_mean_pattern = np.zeros((mesh.nz, mesh.ny, nx_pattern))
        V_mean_pattern = np.zeros((mesh.nz, mesh.ny, nx_pattern))
        W_mean_pattern = np.zeros((mesh.nz, mesh.ny, nx_pattern))
        P_mean_pattern = np.zeros((mesh.nz, mesh.ny, nx_pattern))
        
        if (settings.get_T):
            T_mean_pattern = np.zeros((mesh.nz, mesh.ny, nx_pattern))
        
        for n in range(num_of_pattern):
            
            U_mean_pattern += self.U.mean[:, :, nx_pattern*n:nx_pattern*(n+1)]
            V_mean_pattern += self.V.mean[:, :, nx_pattern*n:nx_pattern*(n+1)]
            W_mean_pattern += self.W.mean[:, :, nx_pattern*n:nx_pattern*(n+1)]
            P_mean_pattern += self.P.mean[:, :, nx_pattern*n:nx_pattern*(n+1)]
            
            if (settings.get_T):
                T_mean_pattern += self.T.mean[:, :, nx_pattern*n:nx_pattern*(n+1)]
    
        U_mean_pattern = U_mean_pattern/num_of_pattern
        V_mean_pattern = V_mean_pattern/num_of_pattern
        W_mean_pattern = W_mean_pattern/num_of_pattern
        P_mean_pattern = P_mean_pattern/num_of_pattern
        
        if (settings.get_T):
            T_mean_pattern = T_mean_pattern/num_of_pattern
            
            T_mean_3D = np.reshape(self.T.mean_spot,(mesh.ny,1))
            T_mean_3D = np.tile(T_mean_3D,(mesh.nz,1,mesh.nx))
    
        U_mean_3D = np.reshape(self.U.mean_spot,(mesh.ny,1))
        U_mean_3D = np.tile(U_mean_3D,(mesh.nz,1,mesh.nx))
    
        V_mean_3D = np.reshape(self.V.mean_spot,(mesh.ny,1))
        V_mean_3D = np.tile(V_mean_3D,(mesh.nz,1,mesh.nx))
    
        W_mean_3D = np.reshape(self.W.mean_spot,(mesh.ny,1))
        W_mean_3D = np.tile(W_mean_3D,(mesh.nz,1,mesh.nx))
    
        P_mean_3D = np.reshape(self.P.mean_spot,(mesh.ny,1))
        P_mean_3D = np.tile(P_mean_3D,(mesh.nz,1,mesh.nx))
        
        for n in range(num_of_pattern):
            
            U_mean[:, :, nx_pattern*n:nx_pattern*(n+1)] = U_mean_pattern - U_mean_3D[:, :, nx_pattern*n:nx_pattern*(n+1)]
            V_mean[:, :, nx_pattern*n:nx_pattern*(n+1)] = V_mean_pattern - V_mean_3D[:, :, nx_pattern*n:nx_pattern*(n+1)]
            W_mean[:, :, nx_pattern*n:nx_pattern*(n+1)] = W_mean_pattern - W_mean_3D[:, :, nx_pattern*n:nx_pattern*(n+1)]
            P_mean[:, :, nx_pattern*n:nx_pattern*(n+1)] = P_mean_pattern - P_mean_3D[:, :, nx_pattern*n:nx_pattern*(n+1)]
            
            if (settings.get_T):
                T_mean[:, :, nx_pattern*n:nx_pattern*(n+1)] = T_mean_pattern - T_mean_3D[:, :, nx_pattern*n:nx_pattern*(n+1)]
    
        del U_mean_3D, V_mean_3D, W_mean_3D, P_mean_3D
        del U_mean_pattern, V_mean_pattern, W_mean_pattern, P_mean_pattern
        
        if (settings.get_T):
            del T_mean_3D
            del T_mean_pattern
    
        self.U.set_mean(U_mean)
        self.V.set_mean(V_mean)
        self.W.set_mean(W_mean)
        self.P.set_mean(P_mean)
        
        if (settings.get_T):
            self.T.set_mean(T_mean)
            
        del U_mean, V_mean, W_mean, P_mean, T_mean
        

    def read_current_fields(self, settings:Settings, current_iteration, folder_path='3D'):
        """Read the instantaneous fields corresponding the current_iteration of the simulation.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        current_iteration : Integer
            Iteration corresponding to the fields read and processed
    
            Set instantaneous fields of the CFD_stats object
        """
        
        results_path = settings.root + '/' + settings.current_path + '/' + settings.current_path + '/Results/' + str(folder_path) + '/field' + str(int(np.rint(current_iteration*1000)))
        
        progress = round(100 * (np.where(settings.all_iterations==current_iteration)[0][0] + 1) / len(settings.all_iterations), 2)
        print(str(progress) + " %")
        
        # U (streamwise velocity here)
        U_instantaneous = read_h5_file(results_path, 'W', settings.streamwise)
        # V
        V_instantaneous = read_h5_file(results_path, 'V', settings.streamwise)
        # U (spanwise velocity here)
        W_instantaneous = read_h5_file(results_path, 'U', settings.streamwise)
        # Pressure
        P_instantaneous = read_h5_file(results_path, 'P', settings.streamwise)
        
        if (settings.get_T):
            self.T.set_instantaneous(read_h5_file(results_path, 'sca1', settings.streamwise))
            
        # from staggered to cell center field
        if (settings.IBM_flag):
            U_instantaneous = interpolate_to_cc_x(U_instantaneous*self.spot)
            V_instantaneous = interpolate_to_cc_y(V_instantaneous*self.spot)
            W_instantaneous = interpolate_to_cc_z(W_instantaneous*self.spot)
        else:
            U_instantaneous = interpolate_to_cc_x(U_instantaneous)
            V_instantaneous = interpolate_to_cc_y(V_instantaneous)
            W_instantaneous = interpolate_to_cc_z(W_instantaneous)
        
        self.U.set_instantaneous(U_instantaneous)
        self.V.set_instantaneous(V_instantaneous)
        self.W.set_instantaneous(W_instantaneous)
        self.P.set_instantaneous(P_instantaneous)
        
        del U_instantaneous, V_instantaneous, W_instantaneous, P_instantaneous
        
    # def from_center_to_bulk(self, settings:Settings):
        
    #     reb = np.mean(self.U.instantaneous[:,:,-1]) * settings.Re
        
    #     self.U.set_instantaneous(self.U.instantaneous * (settings.Re/reb))
    #     self.V.set_instantaneous(self.V.instantaneous * (settings.Re/reb))
    #     self.W.set_instantaneous(self.W.instantaneous * (settings.Re/reb))
        
    #     print("New Reynolds bulk of ", reb)
        
        
    def get_turbulent_spot(self, settings:Settings, mesh:CFD_mesh, spot_sensitivity_study=False, blur_coef=6):
        """From the instantaneous flow field, detect and interface the turbulent spot.
        Only used for transitional flows.
        The spot is then used to compute statistics only in it.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Initialize the spot of the CFD_stats object
        """
        
        ###########################################################
        tmp_dUfdz = np.gradient(self.U.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(self.V.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)

        tmp_dVfdx = np.gradient(self.V.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(self.W.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)
        
        tmp_dUfdy = D1_3Dy(self.U.get_fluctuations(), mesh.Yc)
        tmp_dWfdy = D1_3Dy(self.W.get_fluctuations(), mesh.Yc)
                    
        # Omega_x RMS in the spot
        omega_x_rms_spot = (tmp_dWfdy-tmp_dVfdz)**2
        
        # Omega_y RMS in the spot
        omega_y_rms_spot = (tmp_dUfdz-tmp_dWfdx)**2
        
        # Omega_z RMS in the spot
        omega_z_rms_spot = (tmp_dVfdx-tmp_dUfdy)**2
        
        del tmp_dUfdz, tmp_dVfdz, tmp_dVfdx, tmp_dWfdx, tmp_dUfdy, tmp_dWfdy
        ###########################################################
        
        E_fluctuations = omega_z_rms_spot + omega_y_rms_spot + omega_x_rms_spot
        del omega_x_rms_spot, omega_y_rms_spot, omega_z_rms_spot
        E_fluctuations = gaussian_filter(E_fluctuations, blur_coef)
        self.spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))
                    
        for j in range(mesh.ny):
            E_max = np.max(E_fluctuations[:,j,:])
            
            alpha = 0.05
            self.spot[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
            
            
        
        if (spot_sensitivity_study):
            
            self.spot0d01 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            self.spot0d1 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            self.spot0d2 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
            self.spot_no_filter = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            self.spot_filter3 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
            self.spot_kin_energy = np.zeros((mesh.nz,mesh.ny,mesh.nx))
                        
            for j in range(mesh.ny):
                E_max = np.max(E_fluctuations[:,j,:])            
                    
                alpha = 0.01
                self.spot0d01[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
        
                alpha = 0.1
                self.spot0d1[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
                
                alpha = 0.2
                self.spot0d2[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
                
            E_fluctuations = omega_z_rms_spot + omega_y_rms_spot + omega_x_rms_spot
            
            for j in range(mesh.ny):
                alpha = 0.05
                E_max = np.max(E_fluctuations[:,j,:])
                self.spot_no_filter[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
            
            E_fluctuations = gaussian_filter(E_fluctuations, 3)
            for j in range(mesh.ny):
                alpha = 0.05
                E_max = np.max(E_fluctuations[:,j,:])
                self.spot_filter3[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
                
            
            E_fluctuations = self.U.get_fluctuations()**2 + self.V.get_fluctuations()**2 + self.W.get_fluctuations()**2
            E_fluctuations = gaussian_filter(E_fluctuations, 6)
            for j in range(mesh.ny):
                alpha = 0.05
                E_max = np.max(E_fluctuations[:,j,:])
                self.spot_kin_energy[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
            
        ##################################################################
        # Get a more accurate spot
        ##################################################################

        if (settings.heart_spot or settings.wave_packet_spot):
        
            noisy = self.U.instantaneous[:,0,:]/mesh.Yc[0] - self.U.mean[:,0,:]/mesh.Yc[0]
            
            freq_x_theo = (1/settings.Lx) * np.arange(mesh.nx)
            freq_z_theo = (1/settings.Lz) * np.arange(mesh.nz)

            freq_x_theo_2 = freq_x_theo - np.max(freq_x_theo)/2
            freq_x_theo_2 = np.roll(freq_x_theo_2,-mesh.nx//2)
            
            freq_z_theo_2 = freq_z_theo - np.max(freq_z_theo)/2
            freq_z_theo_2 = np.roll(freq_z_theo_2,-mesh.nz//2)
            
            lambda_x = 2*np.pi/0.9 # filter all scale greater than h
            lambda_z = 2*np.pi/0.5 # filter all scale greater than h
            
            gx = np.exp(-((2*np.pi*freq_x_theo_2)**2*lambda_x**2/24))      # Gauss exp kernel in x direction
            gz = np.exp(-((2*np.pi*freq_z_theo_2)**2*lambda_z**2/24))     # Gauss exp kernel in z direction
            g2d = np.outer(gz, gx)                                  # 2d filter kernel in x-z plane
            
            filtered=np.fft.ifft2((np.fft.fft2(noisy)*g2d)).real                           # 2d filter kernel in x-z plane
                        
            wall_spot = np.zeros((mesh.nz,mesh.nx))
                                                                    
            alpha = 0.6
            wall_spot[:,:] = (filtered[:,:]>=alpha*np.max(filtered))

            if (settings.heart_spot):
            
                for j in range(mesh.ny):
                    self.spot[:,j,:] = wall_spot[:,:]#*self.spot[:,j,:]

            elif (settings.wave_packet_spot):
            
                for j in range(mesh.ny):
                    self.spot[:,j,:] = (1-wall_spot[:,:])*self.spot[:,j,:]
            
    
    def get_spot(self, settings:Settings, mesh:CFD_mesh):
        """From the user settings, define the spot.
        Can be used for all kind of flows. The spot can be the total channel or a part of it.
        The spot is then used to compute statistics only in it.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Initialize the spot of the CFD_stats object
        """
        
        if (settings.total_domain):
            self.spot = np.ones((mesh.nz,mesh.ny,mesh.nx))
            
        else:
            xs = find_nearest(mesh.Xc, settings.x_bounds[0])
            xe = find_nearest(mesh.Xc, settings.x_bounds[1])
            
            ys = find_nearest(mesh.Yc, settings.y_bounds[0])
            ye = find_nearest(mesh.Yc, settings.y_bounds[1])
            
            zs = find_nearest(mesh.Zc, settings.z_bounds[0])
            ze = find_nearest(mesh.Zc, settings.z_bounds[1])
            
            self.spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
            for i in range(xs,xe+1):
                for j in range(ys,ye+1):
                    for k in range(zs,ze+1):
                        self.spot[k,j,i] = 1
                
        
    def get_spot_ibm(self, settings:Settings, mesh:CFD_mesh):
        """From the user settings, define the ibm mask, and thus the spot.
        Can be used for all kind of flows. The spot can be the total channel or a part of it.
        The spot is then used to compute statistics only in it.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Initialize the spot of the CFD_stats object
        """
        
        self.spot = np.ones((mesh.nz,mesh.ny,mesh.nx))
        
        # print("Attention au spot ...")
        # self.spot[:,:,:1272+1] = 0
        # self.spot[:,:,2101:] = 0
        
        for n in range(len(settings.i_start)):
            
            self.spot[settings.k_start[n]:settings.k_end[n]+1,settings.j_start[n]:settings.j_end[n]+1,settings.i_start[n]:settings.i_end[n]+1] = 0
                
    
    def get_turbulent_thermal_spot(self, settings:Settings, mesh:CFD_mesh, blur_coef=6):
        """From the instantaneous flow field, detect and interface the turbulent thermal spot.
        Only used for transitional flows.
        The spot is then used to compute statistics only in it.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Initialize the thermal spot of the CFD_stats object
        """
        
        T_fluctuations_spot = gaussian_filter(self.T.get_fluctuations()**2, blur_coef)

        self.thermal_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
        alpha = 0.05

        for j in range(mesh.ny):
            T_max = np.max(T_fluctuations_spot[:,j,:])
            self.thermal_spot[:,j,:] = (T_fluctuations_spot[:,j,:]>=alpha*T_max)
            
        ##################################################################
        # Get a more accurate spot
        ##################################################################

        if (settings.heart_spot or settings.wave_packet_spot):
        
            noisy = self.T.instantaneous[:,0,:]/mesh.Yc[0] - self.T.mean[:,0,:]/mesh.Yc[0]
            
            freq_x_theo = (1/settings.Lx) * np.arange(mesh.nx)
            freq_z_theo = (1/settings.Lz) * np.arange(mesh.nz)

            freq_x_theo_2 = freq_x_theo - np.max(freq_x_theo)/2
            freq_x_theo_2 = np.roll(freq_x_theo_2,-mesh.nx//2)
            
            freq_z_theo_2 = freq_z_theo - np.max(freq_z_theo)/2
            freq_z_theo_2 = np.roll(freq_z_theo_2,-mesh.nz//2)
            
            lambda_x = 2*np.pi/0.9 # filter all scale greater than h
            lambda_z = 2*np.pi/0.5 # filter all scale greater than h
            
            gx = np.exp(-((2*np.pi*freq_x_theo_2)**2*lambda_x**2/24))      # Gauss exp kernel in x direction
            gz = np.exp(-((2*np.pi*freq_z_theo_2)**2*lambda_z**2/24))     # Gauss exp kernel in z direction
            g2d = np.outer(gz, gx)                                  # 2d filter kernel in x-z plane
            
            filtered=np.fft.ifft2((np.fft.fft2(noisy)*g2d)).real                           # 2d filter kernel in x-z plane
                        
            wall_spot = np.zeros((mesh.nz,mesh.nx))
            
            alpha = 0.6
            print(alpha)
            wall_spot[:,:] = (filtered[:,:]>=alpha*np.max(filtered))

            if (settings.heart_spot):
            
                for j in range(mesh.ny):
                    self.thermal_spot[:,j,:] = wall_spot[:,:]#*self.thermal_spot[:,j,:]

            elif (settings.wave_packet_spot):
            
                for j in range(mesh.ny):
                    self.thermal_spot[:,j,:] = (1-wall_spot[:,:])*self.thermal_spot[:,j,:]
            
    
    def get_thermal_spot(self):
        """Initialize the thermal spot.
        Can be used for all kind of flows. The spot can be the total channel or a part of it.
        The spot is then used to compute statistics only in it.
        """
        
        self.thermal_spot = self.spot
            
            
    def get_wall_shear_stress(self, mesh):
        """Get the fluctuated wall shear stress at the lower wall.
        
        Parameters
        ----------
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
        Returns
        -------
        2D numpy array with shape (nz,nx)
            The fluctuated wall shear stress at the lower wall
        """
        
        # compute inst - mean value
        tmp = (self.U.instantaneous[:,0,:]-self.U.mean[:,0,:])/mesh.Yc[0]
        
        return tmp
            
            
    def get_wall_temperature_flux(self, mesh):
        """Get the fluctuated wall temperature flux at the lower wall.
        
        Parameters
        ----------
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
        Returns
        -------
        2D numpy array with shape (nz,nx)
            The fluctuated wall temperature flux at the lower wall
        """
        
        # compute inst - mean value
        tmp = (self.T.instantaneous[:,0,:]-self.T.mean[:,0,:])/mesh.Yc[0]
        
        return tmp
    
    
    def set_mean_spot(self, settings:Settings, mesh:CFD_mesh):
        """Compute the mean in the spot (1D numpy array with shape (ny)) of all velocities variables.
        Only used for transitional flows.
        Compute also basic statistics numbers (Re_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Set the mean velocities fields in the spot of the CFD_stats object
        """
        
        # Re_tau in the spot
        tmp = self.spot[:,0,:] * self.U.instantaneous[:,0,:]
        tmp[tmp==0] = np.nan
        self.Re_tau_spot = np.sqrt( settings.Re * np.nanmean(tmp) / mesh.Yc[0])
        
        if (settings.symmetry):
            tmp = self.spot[:,-1,:] * self.U.instantaneous[:,-1,:]
            tmp[tmp==0] = np.nan
            self.Re_tau_spot += np.sqrt( settings.Re * np.nanmean(tmp) / mesh.Yc[0])
            self.Re_tau_spot = self.Re_tau_spot/2

        ####################################################
        # U_mean in the spot
        tmp = self.spot * self.U.instantaneous
        tmp[tmp==0] = np.nan
        U_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        # V_mean in the spot
        tmp = self.spot * self.V.instantaneous
        tmp[tmp==0] = np.nan
        V_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        # W_mean in the spot
        tmp = self.spot * self.W.instantaneous
        tmp[tmp==0] = np.nan
        W_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        # P_mean in the spot
        tmp = self.spot * self.P.instantaneous
        tmp[tmp==0] = np.nan
        P_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        self.U.set_mean_spot(U_mean_spot)
        self.V.set_mean_spot(V_mean_spot)
        self.W.set_mean_spot(W_mean_spot)
        self.P.set_mean_spot(P_mean_spot)
        
        del U_mean_spot, V_mean_spot, W_mean_spot, P_mean_spot
        
        
    def set_mean_thermal_spot(self, settings:Settings, mesh:CFD_mesh):
        """Compute the mean in the spot (1D numpy array with shape (ny)) of all temperature variables.
        Only used for transitional flows.
        Compute also basic statistics numbers (T_tau, Nu)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
    
            Set the mean temperature fields in the spot of the CFD_stats object
        """
        
        # T_tau in the spot
        tmp = self.thermal_spot[:,0,:] * self.T.instantaneous[:,0,:]
        tmp[tmp==0] = np.nan
        self.T_tau_spot = abs((np.nanmean(tmp)+settings.deltaT) / (mesh.Yc[0]*self.Re_tau_spot*settings.Pr))
    
        if (settings.symmetry):
            tmp = self.thermal_spot[:,-1,:] * self.T.instantaneous[:,-1,:]
            tmp[tmp==0] = np.nan
            self.T_tau_spot += abs((np.nanmean(tmp)-settings.deltaT) / (mesh.Yc[0]*self.Re_tau_spot*settings.Pr))
            self.T_tau_spot = self.T_tau_spot/2
            
        # Nusselt in the spot
        tmp = self.thermal_spot[:,0,:] * self.T.instantaneous[:,0,:]
        tmp[tmp==0] = np.nan
        self.Nu_spot = abs((np.nanmean(tmp)+settings.deltaT) / (mesh.Yc[0])) * ( 1 / settings.deltaT )
    
        if (settings.symmetry):
            tmp = self.thermal_spot[:,-1,:] * self.T.instantaneous[:,-1,:]
            tmp[tmp==0] = np.nan
            self.Nu_spot += abs((np.nanmean(tmp)-settings.deltaT) / (mesh.Yc[0])) * ( 1 / settings.deltaT )
            self.Nu_spot = self.Nu_spot/2
            
        ####################################################
        # T_mean in the spot
        tmp = self.thermal_spot * self.T.instantaneous
        tmp[tmp==0] = np.nan
        T_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))

        self.T.set_mean_spot(T_mean_spot)
        
        del T_mean_spot
        
        
    def compute_basic_stats(self):
        """Compute basic velocity statistics (1D numpy array with shape (ny)) of the flow in the spot.
        """
        
        # U_fluc*U_fluc in the spot
        self.UfUf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), self.U.get_fluctuations_spot())

        # V_fluc*V_fluc in the spot
        self.VfVf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), self.V.get_fluctuations_spot())

        # W_fluc*W_fluc in the spot
        self.WfWf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_spot(), self.W.get_fluctuations_spot())

        # U_fluc*V_fluc in the spot
        self.UfVf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), self.V.get_fluctuations_spot())

        # U_fluc*W_fluc in the spot
        self.UfWf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), self.W.get_fluctuations_spot())

        # V_fluc*W_fluc in the spot
        self.VfWf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), self.W.get_fluctuations_spot())

        # P_fluc*U_fluc in the spot
        self.PfUf_spot += compute_correlations_spot(self.spot * self.P.get_fluctuations_spot(), self.U.get_fluctuations_spot())

        ####################################################
        # U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), self.U.get_fluctuations_spot(), self.U.get_fluctuations_spot())

        # U_fluc*U_fluc*V_fluc in the spot
        self.UfUfVf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), self.U.get_fluctuations_spot(), self.V.get_fluctuations_spot())

        # U_fluc*V_fluc*V_fluc in the spot
        self.UfVfVf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), self.V.get_fluctuations_spot(), self.V.get_fluctuations_spot())

        # V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), self.V.get_fluctuations_spot(), self.V.get_fluctuations_spot())

        # W_fluc*W_fluc*V_fluc in the spot
        self.WfWfVf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_spot(), self.W.get_fluctuations_spot(), self.V.get_fluctuations_spot())

        # W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_spot(), self.W.get_fluctuations_spot(), self.W.get_fluctuations_spot())

        ####################################################
        # U_fluc*U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUfUf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), self.U.get_fluctuations_spot(), self.U.get_fluctuations_spot(), self.U.get_fluctuations_spot())

        # V_fluc*V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVfVf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), self.V.get_fluctuations_spot(), self.V.get_fluctuations_spot(), self.V.get_fluctuations_spot())

        # W_fluc*W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWfWf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_spot(), self.W.get_fluctuations_spot(), self.W.get_fluctuations_spot(), self.W.get_fluctuations_spot())
        
        
    def compute_mean_streamwise_variations(self, settings:Settings, mesh:CFD_mesh):
        """Compute basic velocity statistics (1D numpy array with shape (ny)) of the flow in the spot.
        """
        
        ################## Re_tau(x)
        self.Re_tau_x = np.sqrt(abs(compute_correlations_x(self.spot * self.U.mean, settings.Re/mesh.Yc[0])))[0,:]

        ################## T_tau(x)
        if ('PSK' in settings.current_path) or ('Ka' in settings.current_path):
            tmp = compute_correlations_x(self.spot * self.T.mean, 1)
            # self.T_tau_x = (tmp[-2,:] - tmp[-1,:]) / ((mesh.Yc[1]-mesh.Yc[0])*self.Re_tau_x*settings.Pr)
            self.T_tau_x = (tmp[-2,:] - tmp[-1,:]) / (mesh.Yc[1]-mesh.Yc[0])
            # self.T_tau_x = 1 * (125/self.Re_tau_x)
            
        else:
            self.T_tau_x = abs(settings.deltaT+compute_correlations_x(self.spot * self.T.mean, 1))[0,:] / (mesh.Yc[0]*self.Re_tau_x*settings.Pr)

        ################## Theta_bulk(x)
        u_theta = np.zeros((mesh.nx))
        u = np.zeros((mesh.nx))
        
        
        if ('PSK' in settings.current_path) or ('Ka' in settings.current_path):
            
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta[:] += self.thermal_spot[k,-1-j,:]*self.U.mean[k,-1-j,:]*self.T.mean[k,-1-j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u[:] += self.thermal_spot[k,-1-j,:]*self.U.mean[k,-1-j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            self.theta_b_x = (u_theta/u)
            
        else:
            
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(self.T.mean[k,j,:]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            self.theta_b_x = (u_theta/u)

        ################## Nu(x)
        u_theta = np.zeros((mesh.nx))
        u = np.zeros((mesh.nx))
        
        
        if ('PSK' in settings.current_path) or ('Ka' in settings.current_path):
            
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta[:] += self.thermal_spot[k,-1-j,:]*self.U.mean[k,-1-j,:]*self.T.mean[k,-1-j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u[:] += self.thermal_spot[k,-1-j,:]*self.U.mean[k,-1-j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            theta_b = (u_theta/u) #* (1/self.T_tau_x)
            
            self.Nu_x = 2*self.Re_tau_x*settings.Pr/theta_b
            
        else:
        
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(self.T.mean[k,j,:]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            theta_b = (u_theta/u) * (1/self.T_tau_x)
            
            self.Nu_x = 2*self.Re_tau_x*settings.Pr/theta_b
            
        ################## Ui(x) et T(x)
        self.U_mean_x = np.mean(self.U.mean, axis=0)
        self.V_mean_x = np.mean(self.V.mean, axis=0)
        self.W_mean_x = np.mean(self.W.mean, axis=0)
        self.T_mean_x = np.mean(self.T.mean, axis=0)
            

    def compute_mean_streamwise_variations_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute basic velocity statistics (1D numpy array with shape (ny)) of the flow in the spot.
        """
        
        ################## Re_tau(x)
        j_e  = settings.j_end[0]
        y_end = settings.y_end[0]
        
        tmp = self.spot[:,0,:] * self.U.mean[:,0,:]/mesh.Yc[0] + (1-self.spot[:,j_e,:]) * self.U.mean[:,j_e+1,:]/(mesh.Yc[j_e+1]-y_end)
        tmp[tmp==0] = np.nan
        correlx = np.nan_to_num(np.nanmean(tmp,axis=(0)))
        
        self.Re_tau_x = np.sqrt( correlx * settings.Re )

        ################## T_tau(x)
        if ('PSK' in settings.current_path) or ('Ka' in settings.current_path):
            tmp = compute_correlations_x(self.spot * self.T.mean, 1)
            self.T_tau_x = (tmp[-2,:] - tmp[-1,:]) / ((mesh.Yc[1]-mesh.Yc[0])*self.Re_tau_x*settings.Pr)
            
        else:
            self.T_tau_x = abs(settings.deltaT+compute_correlations_x(self.spot * self.T.mean, 1))[0,:] / (mesh.Yc[0]*self.Re_tau_x*settings.Pr)

        ################## Nu(x)
        u_theta = np.zeros((mesh.nx))
        u = np.zeros((mesh.nx))
        
        
        if ('PSK' in settings.current_path) or ('Ka' in settings.current_path):
            
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta[:] += self.thermal_spot[k,-1-j,:]*self.U.mean[k,-1-j,:]*self.T.mean[k,-1-j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u[:] += self.thermal_spot[k,-1-j,:]*self.U.mean[k,-1-j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            theta_b = (u_theta/u) #* (1/self.T_tau_x)
            
            self.Nu_x = 2*self.Re_tau_x*settings.Pr/theta_b
            
        else:
        
            for k in range(mesh.nz):
                for j in range(mesh.ny//2):
                    
                    u_theta[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(self.T.mean[k,j,:]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
                    u[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                    
            theta_b = (u_theta/u) * (1/self.T_tau_x)
            
            self.Nu_x = 2*self.Re_tau_x*settings.Pr/theta_b
    
    
    
    def compute_local_stats(self, settings:Settings, mesh:CFD_mesh):
        """Compute basic velocity statistics (1D numpy array with shape (ny)) of the flow in the spot.
        """
        
        # ################## Re_tau(x)
        # Re_tau_x_inst = np.sqrt(abs(compute_correlations_x(self.spot * self.U.instantaneous, settings.Re/mesh.Yc[0])))[0,:]
        # self.Re_tau_x += Re_tau_x_inst

        # ################## T_tau(x)
        # T_tau_x_inst = abs(settings.deltaT+compute_correlations_x(self.spot * self.T.instantaneous, 1))[0,:] / (mesh.Yc[0]*Re_tau_x_inst*settings.Pr)
        # self.T_tau_x += T_tau_x_inst

        # ################## Nu(x)
        # u_theta = np.zeros((mesh.nx))
        # u = np.zeros((mesh.nx))
        
        # for k in range(mesh.nz):
        #     for j in range(mesh.ny//2):
                
        #         u_theta[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(self.T.mean[k,j,:]+settings.deltaT)*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                
        #         u[:] += self.thermal_spot[k,j,:]*self.U.mean[k,j,:]*(mesh.Y[j+1]-mesh.Y[j])*mesh.zstep
                
        # theta_b = (u_theta/u) * (1/T_tau_x_inst)
        
        # self.Nu_x += 2*Re_tau_x_inst*settings.Pr/theta_b

        ################## Ek(x)
        U_fluc = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        V_fluc = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        W_fluc = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        T_fluc = np.zeros((mesh.nz, mesh.ny, mesh.nx))
        
        for i in range(mesh.nx):
            for j in range(mesh.ny):
                U_fluc[:,j,i] = self.U.instantaneous[:,j,i] - self.U_mean_x[j,i]
                V_fluc[:,j,i] = self.V.instantaneous[:,j,i] - self.V_mean_x[j,i]
                W_fluc[:,j,i] = self.W.instantaneous[:,j,i] - self.W_mean_x[j,i]
                T_fluc[:,j,i] = self.T.instantaneous[:,j,i] - self.T_mean_x[j,i]
        
        Ek = 0.5 * ( U_fluc**2 + V_fluc**2 + W_fluc**2 )

        # average along x
        Ek_x_inst = compute_correlations_x(self.spot * Ek, 1)

        ################## TT(x)
        TT = ( T_fluc**2 )

        # average along x
        TT_x_inst = compute_correlations_x(self.spot * TT, 1)

        ################## Omega(x)
        tmp_dUfdz = np.gradient(U_fluc, mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(V_fluc, mesh.Zc, edge_order=2, axis=0)

        tmp_dVfdx = np.gradient(V_fluc, mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(W_fluc, mesh.Xc, edge_order=2, axis=2)
        
        tmp_dUfdy = D1_3Dy(U_fluc, mesh.Yc)
        tmp_dWfdy = D1_3Dy(W_fluc, mesh.Yc)
                    
        # Omega_x
        omega_x_rms_spot = (tmp_dWfdy-tmp_dVfdz)**2
        
        # Omega_y
        omega_y_rms_spot = (tmp_dUfdz-tmp_dWfdx)**2
        
        # Omega_z
        omega_z_rms_spot = (tmp_dVfdx-tmp_dUfdy)**2

        del tmp_dUfdz, tmp_dVfdz, tmp_dVfdx, tmp_dWfdx, tmp_dUfdy, tmp_dWfdy
        Omega = 0.5* (omega_z_rms_spot + omega_y_rms_spot + omega_x_rms_spot)
        del omega_x_rms_spot, omega_y_rms_spot, omega_z_rms_spot

        # average along x
        Omega_x_inst = compute_correlations_x(self.spot * Omega, 1)        
        
        ### compute mean in y
        for j in range(mesh.ny):
            
            self.Ek_x[:] += Ek_x_inst[j,:]*(mesh.Y[j+1]-mesh.Y[j])/2
            self.TT_x[:] += TT_x_inst[j,:]*(mesh.Y[j+1]-mesh.Y[j])/2
            self.Omega_x[:] += Omega_x_inst[j,:]*(mesh.Y[j+1]-mesh.Y[j])/2
        
    
    def compute_basic_stats_triple_decomposition(self):
        """Compute basic velocity statistics (1D numpy array with shape (ny)) of the flow in the spot.
        """
        
        # U_fluc*U_fluc in the spot
        self.UfUf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp())

        # V_fluc*V_fluc in the spot
        self.VfVf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp())

        # W_fluc*W_fluc in the spot
        self.WfWf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp())

        # U_fluc*V_fluc in the spot
        self.UfVf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp())

        # U_fluc*W_fluc in the spot
        self.UfWf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp())

        # V_fluc*W_fluc in the spot
        self.VfWf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp())

        # P_fluc*U_fluc in the spot
        self.PfUf_spot += compute_correlations_spot(self.spot * self.P.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp())

        ####################################################
        # U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp())

        # U_fluc*U_fluc*V_fluc in the spot
        self.UfUfVf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp())

        # U_fluc*V_fluc*V_fluc in the spot
        self.UfVfVf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp())

        # V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp())

        # W_fluc*W_fluc*V_fluc in the spot
        self.WfWfVf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp())

        # W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp())

        ####################################################
        # U_fluc*U_fluc*U_fluc*U_fluc in the spot
        self.UfUfUfUf_spot += compute_correlations_spot(self.spot * self.U.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp(), self.U.get_fluctuations_triple_decomp())

        # V_fluc*V_fluc*V_fluc*V_fluc in the spot
        self.VfVfVfVf_spot += compute_correlations_spot(self.spot * self.V.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp(), self.V.get_fluctuations_triple_decomp())

        # W_fluc*W_fluc*W_fluc*W_fluc in the spot
        self.WfWfWfWf_spot += compute_correlations_spot(self.spot * self.W.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp(), self.W.get_fluctuations_triple_decomp())
        
        
    def compute_basic_thermal_stats(self):
        """Compute basic temperature statistics (1D numpy array with shape (ny)) of the flow in the spot.
        """
        
        ####################################################
        # T_fluc*T_fluc in the spot
        self.TfTf_spot += compute_correlations_spot(self.thermal_spot * self.T.get_fluctuations_spot(), self.T.get_fluctuations_spot())

        # U_fluc*T_fluc in the spot
        self.UfTf_spot += compute_correlations_spot(self.thermal_spot * self.U.get_fluctuations_spot(), self.T.get_fluctuations_spot())

        # V_fluc*T_fluc in the spot
        self.VfTf_spot += compute_correlations_spot(self.thermal_spot * self.V.get_fluctuations_spot(), self.T.get_fluctuations_spot())

        # W_fluc*T_fluc in the spot
        self.WfTf_spot += compute_correlations_spot(self.thermal_spot * self.W.get_fluctuations_spot(), self.T.get_fluctuations_spot())
        
        ####################################################
        # T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTf_spot += compute_correlations_spot(self.thermal_spot * self.T.get_fluctuations_spot(), self.T.get_fluctuations_spot(), self.T.get_fluctuations_spot())
        
        # V_fluc*T_fluc*T_fluc in the spot
        self.VfTfTf_spot += compute_correlations_spot(self.thermal_spot * self.V.get_fluctuations_spot(), self.T.get_fluctuations_spot(), self.T.get_fluctuations_spot())

        ####################################################
        # T_fluc*T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTfTf_spot += compute_correlations_spot(self.thermal_spot * self.T.get_fluctuations_spot(), self.T.get_fluctuations_spot(), self.T.get_fluctuations_spot(), self.T.get_fluctuations_spot())
        
        
    def compute_basic_thermal_stats_triple_decomposition(self):
        """Compute basic temperature statistics (1D numpy array with shape (ny)) of the flow in the spot.
        """
        
        ####################################################
        # T_fluc*T_fluc in the spot
        self.TfTf_spot += compute_correlations_spot(self.thermal_spot * self.T.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp())

        # U_fluc*T_fluc in the spot
        self.UfTf_spot += compute_correlations_spot(self.thermal_spot * self.U.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp())

        # V_fluc*T_fluc in the spot
        self.VfTf_spot += compute_correlations_spot(self.thermal_spot * self.V.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp())

        # W_fluc*T_fluc in the spot
        self.WfTf_spot += compute_correlations_spot(self.thermal_spot * self.W.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp())
        
        ####################################################
        # T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTf_spot += compute_correlations_spot(self.thermal_spot * self.T.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp())
        
        # V_fluc*T_fluc*T_fluc in the spot
        self.VfTfTf_spot += compute_correlations_spot(self.thermal_spot * self.V.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp())

        ####################################################
        # T_fluc*T_fluc*T_fluc*T_fluc in the spot
        self.TfTfTfTf_spot += compute_correlations_spot(self.thermal_spot * self.T.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp(), self.T.get_fluctuations_triple_decomp())
    
        
    def compute_transport_equation_terms(self, settings:Settings, mesh:CFD_mesh):
        """Compute velocities budgets (uu,uv,vv,ww) (1D numpy array with shape (ny)) of the flow in the spot.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        # Compute stats
        ####################################################

        tmp_dUfdy = D1_3Dy(self.U.get_fluctuations_spot(), mesh.Yc)
        tmp_dVfdy = D1_3Dy(self.V.get_fluctuations_spot(), mesh.Yc)
        tmp_dWfdy = D1_3Dy(self.W.get_fluctuations_spot(), mesh.Yc)

        # dUfdy*dUfdy in the spot
        dUfdydUfdy_spot = compute_correlations_spot(self.spot * tmp_dUfdy, tmp_dUfdy)

        # dUfdy*dVfdy in the spot
        dUfdydVfdy_spot = compute_correlations_spot(self.spot * tmp_dUfdy, tmp_dVfdy)

        # dVfdy*dVfdy in the spot
        dVfdydVfdy_spot = compute_correlations_spot(self.spot * tmp_dVfdy, tmp_dVfdy)

        # dWfdy*dWfdy in the spot
        dWfdydWfdy_spot = compute_correlations_spot(self.spot * tmp_dWfdy, tmp_dWfdy)

        del tmp_dUfdy, tmp_dVfdy, tmp_dWfdy

        ####################################################
        tmp_dUfdx = np.gradient(self.U.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        tmp_dVfdx = np.gradient(self.V.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(self.W.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)

        # dUfdx*dUfdx in the spot
        dUfdxdUfdx_spot = compute_correlations_spot(self.spot * tmp_dUfdx, tmp_dUfdx)

        # dUfdx*dVfdx in the spot
        dUfdxdVfdx_spot = compute_correlations_spot(self.spot * tmp_dUfdx, tmp_dVfdx)

        # dVfdx*dVfdx in the spot
        dVfdxdVfdx_spot = compute_correlations_spot(self.spot * tmp_dVfdx, tmp_dVfdx)

        # dWfdx*dWfdx in the spot
        dWfdxdWfdx_spot = compute_correlations_spot(self.spot * tmp_dWfdx, tmp_dWfdx)

        del tmp_dUfdx, tmp_dVfdx, tmp_dWfdx

        ####################################################
        tmp_dUfdz = np.gradient(self.U.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(self.V.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dWfdz = np.gradient(self.W.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)

        # dUfdz*dUfdz in the spot
        dUfdzdUfdz_spot = compute_correlations_spot(self.spot * tmp_dUfdz, tmp_dUfdz)

        # dUfdz*dVfdz in the spot
        dUfdzdVfdz_spot = compute_correlations_spot(self.spot * tmp_dUfdz, tmp_dVfdz)

        # dVfdz*dVfdz in the spot
        dVfdzdVfdz_spot = compute_correlations_spot(self.spot * tmp_dVfdz, tmp_dVfdz)

        # dWfdz*dWfdz in the spot
        dWfdzdWfdz_spot = compute_correlations_spot(self.spot * tmp_dWfdz, tmp_dWfdz)

        del tmp_dUfdz, tmp_dVfdz, tmp_dWfdz
        
        ####################################################
        tmp_dPfdz = np.gradient(self.P.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        
        tmp_dPfdx = np.gradient(self.P.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)

        tmp_dPfdy = D1_3Dy(self.P.get_fluctuations_spot(), mesh.Yc, correct_at_wall=False)


        # Uf*dPfdx in the spot
        UfdPfdx_spot = compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), tmp_dPfdx)

        # Vf*dPfdx in the spot
        VfdPfdx_spot = compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), tmp_dPfdx)

        # Wf*dPfdz in the spot
        WfdPfdz_spot = compute_correlations_spot(self.spot * self.W.get_fluctuations_spot(), tmp_dPfdz)

        # Uf*dPfdy in the spot
        UfdPfdy_spot = compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), tmp_dPfdy)

        # vf*dPfdy in the spot
        VfdPfdy_spot = compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), tmp_dPfdy)

        del tmp_dPfdx, tmp_dPfdz, tmp_dPfdy

        #dUfUf/dt
        self.uu_budget['Pressure velocity gradient correlation'] += - 2 * UfdPfdx_spot
        self.uu_budget['Dissipation'] += - 2/settings.Re * (dUfdxdUfdx_spot + dUfdydUfdy_spot + dUfdzdUfdz_spot)

        #dUfVf/dt
        self.uv_budget['Pressure velocity gradient correlation'] += - UfdPfdy_spot - VfdPfdx_spot
        self.uv_budget['Dissipation'] += - 2/settings.Re * (dUfdxdVfdx_spot + dUfdydVfdy_spot + dUfdzdVfdz_spot)

        #dVfVf/dt
        self.vv_budget['Pressure velocity gradient correlation'] += - 2 * VfdPfdy_spot
        self.vv_budget['Dissipation'] += - 2/settings.Re * (dVfdxdVfdx_spot + dVfdydVfdy_spot + dVfdzdVfdz_spot)

        #dWfWf/dt
        self.ww_budget['Pressure velocity gradient correlation'] += - 2 * WfdPfdz_spot
        self.ww_budget['Dissipation'] += - 2/settings.Re * (dWfdxdWfdx_spot + dWfdydWfdy_spot + dWfdzdWfdz_spot)
    
     
    def compute_transport_equation_terms_sedat(self, settings:Settings, mesh:CFD_mesh):
        """Compute velocities budgets (uu,uv,vv,ww) (1D numpy array with shape (ny)) of the flow in the spot.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        # Compute stats
        ####################################################

        tmp_dUfdy = D1_3Dy(self.spot * self.U.get_fluctuations_spot(), mesh.Yc)
        tmp_dVfdy = D1_3Dy(self.spot * self.V.get_fluctuations_spot(), mesh.Yc)
        tmp_dWfdy = D1_3Dy(self.spot * self.W.get_fluctuations_spot(), mesh.Yc)
        
        tmp_dUfdy = D1_3Dy_IBM(tmp_dUfdy, self.U.get_fluctuations_spot(), settings, mesh.Yc)
        tmp_dVfdy = D1_3Dy_IBM(tmp_dVfdy, self.V.get_fluctuations_spot(), settings, mesh.Yc)
        tmp_dWfdy = D1_3Dy_IBM(tmp_dWfdy, self.W.get_fluctuations_spot(), settings, mesh.Yc)

        # dUfdy*dUfdy in the spot
        dUfdydUfdy_spot = compute_correlations_spot(self.spot * tmp_dUfdy, tmp_dUfdy)

        # dUfdy*dVfdy in the spot
        dUfdydVfdy_spot = compute_correlations_spot(self.spot * tmp_dUfdy, tmp_dVfdy)

        # dVfdy*dVfdy in the spot
        dVfdydVfdy_spot = compute_correlations_spot(self.spot * tmp_dVfdy, tmp_dVfdy)

        # dWfdy*dWfdy in the spot
        dWfdydWfdy_spot = compute_correlations_spot(self.spot * tmp_dWfdy, tmp_dWfdy)

        del tmp_dUfdy, tmp_dVfdy, tmp_dWfdy

        ####################################################
        tmp_dUfdx = np.gradient(self.spot * self.U.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        tmp_dVfdx = np.gradient(self.spot * self.V.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(self.spot * self.W.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)

        # dUfdx*dUfdx in the spot
        dUfdxdUfdx_spot = compute_correlations_spot(self.spot * tmp_dUfdx, tmp_dUfdx)

        # dUfdx*dVfdx in the spot
        dUfdxdVfdx_spot = compute_correlations_spot(self.spot * tmp_dUfdx, tmp_dVfdx)

        # dVfdx*dVfdx in the spot
        dVfdxdVfdx_spot = compute_correlations_spot(self.spot * tmp_dVfdx, tmp_dVfdx)

        # dWfdx*dWfdx in the spot
        dWfdxdWfdx_spot = compute_correlations_spot(self.spot * tmp_dWfdx, tmp_dWfdx)

        del tmp_dUfdx, tmp_dVfdx, tmp_dWfdx

        ####################################################
        tmp_dUfdz = np.gradient(self.spot * self.U.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(self.spot * self.V.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dWfdz = np.gradient(self.spot * self.W.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)

        # dUfdz*dUfdz in the spot
        dUfdzdUfdz_spot = compute_correlations_spot(self.spot * tmp_dUfdz, tmp_dUfdz)

        # dUfdz*dVfdz in the spot
        dUfdzdVfdz_spot = compute_correlations_spot(self.spot * tmp_dUfdz, tmp_dVfdz)

        # dVfdz*dVfdz in the spot
        dVfdzdVfdz_spot = compute_correlations_spot(self.spot * tmp_dVfdz, tmp_dVfdz)

        # dWfdz*dWfdz in the spot
        dWfdzdWfdz_spot = compute_correlations_spot(self.spot * tmp_dWfdz, tmp_dWfdz)

        del tmp_dUfdz, tmp_dVfdz, tmp_dWfdz
        
        ####################################################
        tmp_dPfdz = np.gradient(self.spot * self.P.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        
        tmp_dPfdx = np.gradient(self.spot * self.P.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)

        tmp_dPfdy = D1_3Dy(self.spot * self.P.get_fluctuations_spot(), mesh.Yc, correct_at_wall=False)
        
        tmp_dPfdy = D1_3Dy_IBM(tmp_dPfdy, self.P.get_fluctuations_spot(), settings, mesh.Yc)


        # Uf*dPfdx in the spot
        UfdPfdx_spot = compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), tmp_dPfdx)

        # Vf*dPfdx in the spot
        VfdPfdx_spot = compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), tmp_dPfdx)

        # Wf*dPfdz in the spot
        WfdPfdz_spot = compute_correlations_spot(self.spot * self.W.get_fluctuations_spot(), tmp_dPfdz)

        # Uf*dPfdy in the spot
        UfdPfdy_spot = compute_correlations_spot(self.spot * self.U.get_fluctuations_spot(), tmp_dPfdy)

        # vf*dPfdy in the spot
        VfdPfdy_spot = compute_correlations_spot(self.spot * self.V.get_fluctuations_spot(), tmp_dPfdy)

        del tmp_dPfdx, tmp_dPfdz, tmp_dPfdy

        #dUfUf/dt
        self.uu_budget['Pressure velocity gradient correlation'] += - 2 * UfdPfdx_spot
        self.uu_budget['Dissipation'] += - 2/settings.Re * (dUfdxdUfdx_spot + dUfdydUfdy_spot + dUfdzdUfdz_spot)

        #dUfVf/dt
        self.uv_budget['Pressure velocity gradient correlation'] += - UfdPfdy_spot - VfdPfdx_spot
        self.uv_budget['Dissipation'] += - 2/settings.Re * (dUfdxdVfdx_spot + dUfdydVfdy_spot + dUfdzdVfdz_spot)

        #dVfVf/dt
        self.vv_budget['Pressure velocity gradient correlation'] += - 2 * VfdPfdy_spot
        self.vv_budget['Dissipation'] += - 2/settings.Re * (dVfdxdVfdx_spot + dVfdydVfdy_spot + dVfdzdVfdz_spot)

        #dWfWf/dt
        self.ww_budget['Pressure velocity gradient correlation'] += - 2 * WfdPfdz_spot
        self.ww_budget['Dissipation'] += - 2/settings.Re * (dWfdxdWfdx_spot + dWfdydWfdy_spot + dWfdzdWfdz_spot)    
        
        
    def compute_transport_equation_terms_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute velocities budgets (uu,uv,vv,ww) (1D numpy array with shape (ny)) of the flow in the spot.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        # Compute stats
        ####################################################

        tmp_dUfdy = D1_3Dy(self.U.get_fluctuations(), mesh.Yc)
        tmp_dVfdy = D1_3Dy(self.V.get_fluctuations(), mesh.Yc)
        tmp_dWfdy = D1_3Dy(self.W.get_fluctuations(), mesh.Yc)

        # dUfdy*dUfdy in the spot
        dUfdydUfdy_spot = compute_correlations_spot(self.spot * tmp_dUfdy, tmp_dUfdy)

        # dUfdy*dVfdy in the spot
        dUfdydVfdy_spot = compute_correlations_spot(self.spot * tmp_dUfdy, tmp_dVfdy)

        # dVfdy*dVfdy in the spot
        dVfdydVfdy_spot = compute_correlations_spot(self.spot * tmp_dVfdy, tmp_dVfdy)

        # dWfdy*dWfdy in the spot
        dWfdydWfdy_spot = compute_correlations_spot(self.spot * tmp_dWfdy, tmp_dWfdy)

        del tmp_dUfdy, tmp_dVfdy, tmp_dWfdy

        ####################################################
        tmp_dUfdx = np.gradient(self.U.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)
        tmp_dVfdx = np.gradient(self.V.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(self.W.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)

        # dUfdx*dUfdx in the spot
        dUfdxdUfdx_spot = compute_correlations_spot(self.spot * tmp_dUfdx, tmp_dUfdx)

        # dUfdx*dVfdx in the spot
        dUfdxdVfdx_spot = compute_correlations_spot(self.spot * tmp_dUfdx, tmp_dVfdx)

        # dVfdx*dVfdx in the spot
        dVfdxdVfdx_spot = compute_correlations_spot(self.spot * tmp_dVfdx, tmp_dVfdx)

        # dWfdx*dWfdx in the spot
        dWfdxdWfdx_spot = compute_correlations_spot(self.spot * tmp_dWfdx, tmp_dWfdx)

        del tmp_dUfdx, tmp_dVfdx, tmp_dWfdx

        ####################################################
        tmp_dUfdz = np.gradient(self.U.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(self.V.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)
        tmp_dWfdz = np.gradient(self.W.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)

        # dUfdz*dUfdz in the spot
        dUfdzdUfdz_spot = compute_correlations_spot(self.spot * tmp_dUfdz, tmp_dUfdz)

        # dUfdz*dVfdz in the spot
        dUfdzdVfdz_spot = compute_correlations_spot(self.spot * tmp_dUfdz, tmp_dVfdz)

        # dVfdz*dVfdz in the spot
        dVfdzdVfdz_spot = compute_correlations_spot(self.spot * tmp_dVfdz, tmp_dVfdz)

        # dWfdz*dWfdz in the spot
        dWfdzdWfdz_spot = compute_correlations_spot(self.spot * tmp_dWfdz, tmp_dWfdz)

        del tmp_dUfdz, tmp_dVfdz, tmp_dWfdz
        
        ####################################################
        tmp_dPfdz = np.gradient(self.P.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)
        
        tmp_dPfdx = np.gradient(self.P.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)

        tmp_dPfdy = D1_3Dy(self.P.get_fluctuations(), mesh.Yc, correct_at_wall=False)


        # Uf*dPfdx in the spot
        UfdPfdx_spot = compute_correlations_spot(self.spot * self.U.get_fluctuations(), tmp_dPfdx)

        # Vf*dPfdx in the spot
        VfdPfdx_spot = compute_correlations_spot(self.spot * self.V.get_fluctuations(), tmp_dPfdx)

        # Wf*dPfdz in the spot
        WfdPfdz_spot = compute_correlations_spot(self.spot * self.W.get_fluctuations(), tmp_dPfdz)

        # Uf*dPfdy in the spot
        UfdPfdy_spot = compute_correlations_spot(self.spot * self.U.get_fluctuations(), tmp_dPfdy)

        # vf*dPfdy in the spot
        VfdPfdy_spot = compute_correlations_spot(self.spot * self.V.get_fluctuations(), tmp_dPfdy)

        del tmp_dPfdx, tmp_dPfdz, tmp_dPfdy

        #dUfUf/dt
        self.uu_budget['Pressure velocity gradient correlation'] += - 2 * UfdPfdx_spot
        self.uu_budget['Dissipation'] += - 2/settings.Re * (dUfdxdUfdx_spot + dUfdydUfdy_spot + dUfdzdUfdz_spot)

        #dUfVf/dt
        self.uv_budget['Pressure velocity gradient correlation'] += - UfdPfdy_spot - VfdPfdx_spot
        self.uv_budget['Dissipation'] += - 2/settings.Re * (dUfdxdVfdx_spot + dUfdydVfdy_spot + dUfdzdVfdz_spot)

        #dVfVf/dt
        self.vv_budget['Pressure velocity gradient correlation'] += - 2 * VfdPfdy_spot
        self.vv_budget['Dissipation'] += - 2/settings.Re * (dVfdxdVfdx_spot + dVfdydVfdy_spot + dVfdzdVfdz_spot)

        #dWfWf/dt
        self.ww_budget['Pressure velocity gradient correlation'] += - 2 * WfdPfdz_spot
        self.ww_budget['Dissipation'] += - 2/settings.Re * (dWfdxdWfdx_spot + dWfdydWfdy_spot + dWfdzdWfdz_spot)
        
    def compute_premultiplied_spectrums(self, settings:Settings, mesh:CFD_mesh, TC=False):
        """Compute premultiplied spectrums (uu,uv,vv,ww) (3D numpy array with shape (nz,ny/2,nx)) of the flow in the spot.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        if (TC):
                
            # Only use in TC of the spot
            filtre_de_hann_z = np.zeros((mesh.nz,mesh.nx))
            filtre_de_hann_x = np.zeros((mesh.nz,mesh.nx))
            
            shift_index = settings.shift

            spot_centerline = np.roll(self.spot[:,0,:],-shift_index,axis=-1)        

            width_x = np.sum(spot_centerline,axis=1)
            width_z = np.sum(spot_centerline,axis=0)

            max_width_x = np.max(width_x)
            max_width_z = np.max(width_z)
            
            i_start = np.where(spot_centerline[mesh.nz//2,:]==1)[0][0]
            i_end = np.where(spot_centerline[mesh.nz//2,:]==1)[0][-1]
            
            
            k_start = np.where(spot_centerline[:,i_start]==1)[0][0]
            k_end = np.where(spot_centerline[:,i_start]==1)[0][-1]

            for i in range(i_start,i_end+1):
                
                filtre_de_hann_x[:,i] = 0.5 - 0.5 * np.cos(2*np.pi*(i-i_start)/(max_width_x))

            for k in range(k_start,k_end+1):
                
                filtre_de_hann_z[k,:] = 0.5 - 0.5 * np.cos(2*np.pi*(k-k_start)/max_width_z)
                
            window_function = filtre_de_hann_x * filtre_de_hann_z * spot_centerline
            window_function = np.roll(window_function,shift_index,axis=-1)
            
            self.window_function = window_function
            
            self.length_x = np.arange(1,np.floor(int(max_width_x)/2),dtype="int")
            mesh.width_x = max_width_x
        
            # Then compute spectrums ...
            self.Lx_TC = mesh.Xc[int(max_width_x)]
            self.Lz_TC = mesh.Zc[int(max_width_z)]
            
            self.width_x = max_width_x
            self.width_z = max_width_z
            
            # The size of the array for the spectrums has to be redefined
            self.uu_spectrum = np.zeros((int(max_width_z)//2-1,mesh.ny//2,int(max_width_x)//2-1))
            self.vv_spectrum = np.zeros((int(max_width_z)//2-1,mesh.ny//2,int(max_width_x)//2-1))
            self.ww_spectrum = np.zeros((int(max_width_z)//2-1,mesh.ny//2,int(max_width_x)//2-1))
            self.uv_spectrum = np.zeros((int(max_width_z)//2-1,mesh.ny//2,int(max_width_x)//2-1))
            
    
            for j in range(mesh.ny//2):
            
                # u fluctuations
                tmp_slice = np.roll(self.U.get_fluctuations_spot()[:,j,:]*self.window_function,-shift_index,axis=-1)[k_start:k_end+1,i_start:i_end+1]
                fhat_2D_u = np.fft.fft2(tmp_slice)[:,self.length_x]
                fhat_2D_u = fhat_2D_u / (int(max_width_x)*int(max_width_z))
                
                # v fluctuations
                tmp_slice = np.roll(self.V.get_fluctuations_spot()[:,j,:]*self.window_function,-shift_index,axis=-1)[k_start:k_end+1,i_start:i_end+1]
                fhat_2D_v = np.fft.fft2(tmp_slice)[:,self.length_x]
                fhat_2D_v = fhat_2D_v / (int(max_width_x)*int(max_width_z))
                
                # w fluctuations
                tmp_slice = np.roll(self.W.get_fluctuations_spot()[:,j,:]*self.window_function,-shift_index,axis=-1)[k_start:k_end+1,i_start:i_end+1]
                fhat_2D_w = np.fft.fft2(tmp_slice)[:,self.length_x]
                fhat_2D_w = fhat_2D_w / (int(max_width_x)*int(max_width_z))
            
                for i in range(int(max_width_x)//2-1):
                    for k in range(int(max_width_z)//2-1):
                        
                        # premultiplied uu spectrum
                        tmp1 = fhat_2D_u[k,i]
                        tmp2 = fhat_2D_u[int(max_width_z)-k-1,i]
                        
                        self.uu_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_TC) * (2*np.pi*i/self.Lx_TC)
                            
                        # premultiplied vv spectrum
                        tmp1 = fhat_2D_v[k,i]
                        tmp2 = fhat_2D_v[int(max_width_z)-k-1,i]
                        
                        self.vv_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_TC) * (2*np.pi*i/self.Lx_TC)
                            
                        # premultiplied ww spectrum
                        tmp1 = fhat_2D_w[k,i]
                        tmp2 = fhat_2D_w[int(max_width_z)-k-1,i]
                        
                        self.ww_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_TC) * (2*np.pi*i/self.Lx_TC)
                            
                        # premultiplied ww spectrum
                        tmp1_u = fhat_2D_u[k,i]
                        tmp2_u = fhat_2D_u[int(max_width_z)-k-1,i]
                        
                        tmp1_v = fhat_2D_v[k,i]
                        tmp2_v = fhat_2D_v[int(max_width_z)-k-1,i]
                        
                        self.uv_spectrum[k,j,i] += ( abs(tmp1_v * np.conj(tmp1_u)) + abs(tmp2_v * np.conj(tmp2_u)) )* (2*np.pi*k/self.Lz_TC) * (2*np.pi*i/self.Lx_TC)
    
            del fhat_2D_u, fhat_2D_v, fhat_2D_w
            
            
        else:
        
            for j in range(mesh.ny//2):
            
                # u fluctuations
                fhat_2D_u = np.fft.fft2((self.U.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_u = fhat_2D_u / (mesh.nx*mesh.nz)
                
                # v fluctuations
                fhat_2D_v = np.fft.fft2((self.V.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_v = fhat_2D_v / (mesh.nx*mesh.nz)
                
                # w fluctuations
                fhat_2D_w = np.fft.fft2((self.W.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_w = fhat_2D_w / (mesh.nx*mesh.nz)
            
                for i in range(mesh.nx//2-1):
                    for k in range(mesh.nz//2-1):
                        
                        # premultiplied uu spectrum
                        tmp1 = fhat_2D_u[k,i]
                        tmp2 = fhat_2D_u[mesh.nz-k-1,i]
                        
                        self.uu_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
                            
                        # premultiplied vv spectrum
                        tmp1 = fhat_2D_v[k,i]
                        tmp2 = fhat_2D_v[mesh.nz-k-1,i]
                        
                        self.vv_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
                            
                        # premultiplied ww spectrum
                        tmp1 = fhat_2D_w[k,i]
                        tmp2 = fhat_2D_w[mesh.nz-k-1,i]
                        
                        self.ww_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
                            
                        # premultiplied ww spectrum
                        tmp1_u = fhat_2D_u[k,i]
                        tmp2_u = fhat_2D_u[mesh.nz-k-1,i]
                        
                        tmp1_v = fhat_2D_v[k,i]
                        tmp2_v = fhat_2D_v[mesh.nz-k-1,i]
                        
                        self.uv_spectrum[k,j,i] += ( abs(tmp1_v * np.conj(tmp1_u)) + abs(tmp2_v * np.conj(tmp2_u)) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
    
            del fhat_2D_u, fhat_2D_v, fhat_2D_w
            
            
    def init_premultiplied_spectrums_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute premultiplied spectrums (uu,uv,vv,ww) (3D numpy array with shape (nz,ny/2,nx)) of the flow while using the immersed boundary method in the domain without the object.
        The idea is apply Hann windowing function for 0<y*<k* on half-pattern and classic FFT computation for k*<y*<1
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        ### First, define the length and width of the pattern
        self.width_x = 0
        self.width_z = 0
        self.number_of_patterns = int(len(settings.x_start)/4 - 1) # Caution, in reality there is one more pattern, it was removed for sake of simplicity

        for n in range(self.number_of_patterns):

            # First half of the channel in the spanwise direction
            i_start_tmp = settings.i_end[int(4*n)] + 1 
            i_end_tmp = settings.i_start[int(4*(n+1))] - 1

            self.width_x = max(self.width_x, i_end_tmp-i_start_tmp) # Caution, this depends on how the roughness are placed ...
            self.width_z = mesh.nz//2 - 1

        self.width_y = settings.j_end[0]

        self.Lx_pattern = mesh.Xc[int(self.width_x)]
        self.Lz_pattern = mesh.Zc[int(self.width_z)]
        self.Ly_pattern = settings.y_end[0]

        ### Then compute the hann filter
        filtre_de_hann_z = np.zeros((self.width_z,self.width_x))
        filtre_de_hann_x = np.zeros((self.width_z,self.width_x))

        i_start_tmp = settings.i_end[0] + 1
        i_end_tmp = settings.i_start[4] - 1

        # reset the start at 0
        i_end_tmp = i_end_tmp - i_start_tmp
        i_start_tmp = i_start_tmp - i_start_tmp

        k_start_tmp = 0
        k_end_tmp = self.width_z

        for i in range(i_start_tmp,i_end_tmp):
            
            filtre_de_hann_x[:,i] = 0.5 - 0.5 * np.cos(2*np.pi*(i-i_start_tmp)/(self.width_x))

        for k in range(k_start_tmp,k_end_tmp):
            
            filtre_de_hann_z[k,:] = 0.5 - 0.5 * np.cos(2*np.pi*(k-k_start_tmp)/self.width_z)
                
        self.window_function = filtre_de_hann_x * filtre_de_hann_z
            
        self.length_x_pattern = np.arange(1,np.floor(int(self.width_x)/2),dtype="int")
                        
        ### The size of the array for the spectrums has to be redefined
        # for 0<y*<k*
        self.uu_spectrum_pattern = np.zeros((int(self.width_z)//2-1,self.width_y,int(self.width_x)//2-1))
        self.vv_spectrum_pattern = np.zeros((int(self.width_z)//2-1,self.width_y,int(self.width_x)//2-1))
        self.ww_spectrum_pattern = np.zeros((int(self.width_z)//2-1,self.width_y,int(self.width_x)//2-1))
        self.uv_spectrum_pattern = np.zeros((int(self.width_z)//2-1,self.width_y,int(self.width_x)//2-1))
        
        if (settings.get_T):
            self.tt_spectrum_pattern = np.zeros((int(self.width_z)//2-1,self.width_y,int(self.width_x)//2-1))


    def compute_premultiplied_spectrums_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute premultiplied spectrums (uu,uv,vv,ww) (3D numpy array with shape (nz,ny/2,nx)) of the flow while using the immersed boundary method in the domain without the object.
        The idea is apply Hann windowing function for 0<y*<k* on half-pattern and classic FFT computation for k*<y*<1
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        ### Finally, compute the spectrums
        for j in range(mesh.ny//2):

            if (j<self.width_y): # for 0<y*<k*

                for n in range(self.number_of_patterns):

                    # First half of the channel in the spanwise direction
                    i_start_tmp = settings.i_end[int(4*n)] + 1
                    i_end_tmp = i_start_tmp + self.width_x

                    k_start_tmp = settings.k_start[int(4*n)]
                    k_end_tmp = k_start_tmp + self.width_z

                    # u fluctuations
                    fhat_2D_u = np.fft.fft2(self.U.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_u = fhat_2D_u / (int(self.width_x)*int(self.width_z))
                    
                    # v fluctuations
                    fhat_2D_v = np.fft.fft2(self.V.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_v = fhat_2D_v / (int(self.width_x)*int(self.width_z))
                    
                    # w fluctuations
                    fhat_2D_w = np.fft.fft2(self.W.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_w = fhat_2D_w / (int(self.width_x)*int(self.width_z))
                
                    for i in range(int(self.width_x)//2-1):
                        for k in range(int(self.width_z)//2-1):
                            
                            # premultiplied uu spectrum
                            tmp1 = fhat_2D_u[k,i]
                            tmp2 = fhat_2D_u[int(self.width_z)-k-1,i]
                            
                            self.uu_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)
                                
                            # premultiplied vv spectrum
                            tmp1 = fhat_2D_v[k,i]
                            tmp2 = fhat_2D_v[int(self.width_z)-k-1,i]
                            
                            self.vv_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)
                                
                            # premultiplied ww spectrum
                            tmp1 = fhat_2D_w[k,i]
                            tmp2 = fhat_2D_w[int(self.width_z)-k-1,i]
                            
                            self.ww_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)
                                
                            # premultiplied ww spectrum
                            tmp1_u = fhat_2D_u[k,i]
                            tmp2_u = fhat_2D_u[int(self.width_z)-k-1,i]
                            
                            tmp1_v = fhat_2D_v[k,i]
                            tmp2_v = fhat_2D_v[int(self.width_z)-k-1,i]
                            
                            self.uv_spectrum_pattern[k,j,i] += ( abs(tmp1_v * np.conj(tmp1_u)) + abs(tmp2_v * np.conj(tmp2_u)) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)


                    # Second half of the channel in the spanwise direction
                    i_start_tmp = settings.i_end[int(4*n+1)] + 1
                    i_end_tmp = i_start_tmp + self.width_x

                    k_start_tmp = settings.k_start[int(4*n+1)]
                    k_end_tmp = k_start_tmp + self.width_z

                    # u fluctuations
                    fhat_2D_u = np.fft.fft2(self.U.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_u = fhat_2D_u / (int(self.width_x)*int(self.width_z))
                    
                    # v fluctuations
                    fhat_2D_v = np.fft.fft2(self.V.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_v = fhat_2D_v / (int(self.width_x)*int(self.width_z))
                    
                    # w fluctuations
                    fhat_2D_w = np.fft.fft2(self.W.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_w = fhat_2D_w / (int(self.width_x)*int(self.width_z))
                
                    for i in range(int(self.width_x)//2-1):
                        for k in range(int(self.width_z)//2-1):
                            
                            # premultiplied uu spectrum
                            tmp1 = fhat_2D_u[k,i]
                            tmp2 = fhat_2D_u[int(self.width_z)-k-1,i]
                            
                            self.uu_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)
                                
                            # premultiplied vv spectrum
                            tmp1 = fhat_2D_v[k,i]
                            tmp2 = fhat_2D_v[int(self.width_z)-k-1,i]
                            
                            self.vv_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)
                                
                            # premultiplied ww spectrum
                            tmp1 = fhat_2D_w[k,i]
                            tmp2 = fhat_2D_w[int(self.width_z)-k-1,i]
                            
                            self.ww_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)
                                
                            # premultiplied ww spectrum
                            tmp1_u = fhat_2D_u[k,i]
                            tmp2_u = fhat_2D_u[int(self.width_z)-k-1,i]
                            
                            tmp1_v = fhat_2D_v[k,i]
                            tmp2_v = fhat_2D_v[int(self.width_z)-k-1,i]
                            
                            self.uv_spectrum_pattern[k,j,i] += ( abs(tmp1_v * np.conj(tmp1_u)) + abs(tmp2_v * np.conj(tmp2_u)) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)

            else: # for k*<y*<1

                # u fluctuations
                fhat_2D_u = np.fft.fft2((self.U.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_u = fhat_2D_u / (mesh.nx*mesh.nz)
                
                # v fluctuations
                fhat_2D_v = np.fft.fft2((self.V.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_v = fhat_2D_v / (mesh.nx*mesh.nz)
                
                # w fluctuations
                fhat_2D_w = np.fft.fft2((self.W.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_w = fhat_2D_w / (mesh.nx*mesh.nz)
            
                for i in range(mesh.nx//2-1):
                    for k in range(mesh.nz//2-1):
                        
                        # premultiplied uu spectrum
                        tmp1 = fhat_2D_u[k,i]
                        tmp2 = fhat_2D_u[mesh.nz-k-1,i]
                        
                        self.uu_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
                            
                        # premultiplied vv spectrum
                        tmp1 = fhat_2D_v[k,i]
                        tmp2 = fhat_2D_v[mesh.nz-k-1,i]
                        
                        self.vv_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
                            
                        # premultiplied ww spectrum
                        tmp1 = fhat_2D_w[k,i]
                        tmp2 = fhat_2D_w[mesh.nz-k-1,i]
                        
                        self.ww_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
                            
                        # premultiplied ww spectrum
                        tmp1_u = fhat_2D_u[k,i]
                        tmp2_u = fhat_2D_u[mesh.nz-k-1,i]
                        
                        tmp1_v = fhat_2D_v[k,i]
                        tmp2_v = fhat_2D_v[mesh.nz-k-1,i]
                        
                        self.uv_spectrum[k,j,i] += ( abs(tmp1_v * np.conj(tmp1_u)) + abs(tmp2_v * np.conj(tmp2_u)) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)


        del fhat_2D_u, fhat_2D_v, fhat_2D_w        
        
        
    def save_transport_equation_terms(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save velocities budgets (uu,uv,vv,ww) to file.
        Compute reamining quantities needed for the budgets.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        ####################################################
        #Computation of all the gradients needed

        #gradient mean_velocity
        dUmeandy = D1_3Dy(self.U.mean_spot, mesh.Yc)

        dUfUfVf_Meandy = D1_3Dy(self.UfUfVf_spot, mesh.Yc)
        dUfVfVf_Meandy = D1_3Dy(self.UfVfVf_spot, mesh.Yc)
        dVfVfVf_Meandy = D1_3Dy(self.VfVfVf_spot, mesh.Yc)
        dWfWfVf_Meandy = D1_3Dy(self.WfWfVf_spot, mesh.Yc)

        d2UfUf_Meandy2 = D2_3Dy(self.UfUf_spot, mesh.Yc, mesh.Y)
        d2UfVf_Meandy2 = D2_3Dy(self.UfVf_spot, mesh.Yc, mesh.Y)
        d2VfVf_Meandy2 = D2_3Dy(self.VfVf_spot, mesh.Yc, mesh.Y)
        d2WfWf_Meandy2 = D2_3Dy(self.WfWf_spot, mesh.Yc, mesh.Y)
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TransportEquationsTerms/'
        createPostProcessingDirectory(current_path_postProcess)

        #dUfUf/dt
        self.uu_budget['Production'] = - 2*self.UfVf_spot * dUmeandy
        self.uu_budget['Turbulent Diffusive Flux'] = - dUfUfVf_Meandy
        self.uu_budget['Molecular diffusion'] = 1/settings.Re * d2UfUf_Meandy2
                
        write_stats(self.uu_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dUfVf/dt
        self.uv_budget['Production'] = - self.VfVf_spot * dUmeandy
        self.uv_budget['Turbulent Diffusive Flux'] = - dUfVfVf_Meandy
        self.uv_budget['Molecular diffusion'] = 1/settings.Re * d2UfVf_Meandy2
                
        write_stats(self.uv_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dVfVf/dt
        self.vv_budget['Turbulent Diffusive Flux'] = - dVfVfVf_Meandy
        self.vv_budget['Molecular diffusion'] = 1/settings.Re * d2VfVf_Meandy2
                
        write_stats(self.vv_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dWfWf/dt
        self.ww_budget['Turbulent Diffusive Flux'] = - dWfWfVf_Meandy
        self.ww_budget['Molecular diffusion'] = 1/settings.Re * d2WfWf_Meandy2
                
        write_stats(self.ww_budget, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def save_transport_equation_terms_sedat(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save velocities budgets (uu,uv,vv,ww) to file.
        Compute reamining quantities needed for the budgets.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        ####################################################
        #Computation of all the gradients needed

        #gradient mean_velocity
        
        U_mean_3D = np.reshape(self.U.mean_spot,(mesh.ny,1))
        U_mean_3D = np.tile(U_mean_3D,(mesh.nz,1,mesh.nx))
        
        #
        UfUfVf_3D = np.reshape(self.UfUfVf_spot,(mesh.ny,1))
        UfUfVf_3D = np.tile(UfUfVf_3D,(mesh.nz,1,mesh.nx))
        
        UfVfVf_3D = np.reshape(self.UfVfVf_spot,(mesh.ny,1))
        UfVfVf_3D = np.tile(UfVfVf_3D,(mesh.nz,1,mesh.nx))
        
        VfVfVf_3D = np.reshape(self.VfVfVf_spot,(mesh.ny,1))
        VfVfVf_3D = np.tile(VfVfVf_3D,(mesh.nz,1,mesh.nx))
        
        WfWfVf_3D = np.reshape(self.WfWfVf_spot,(mesh.ny,1))
        WfWfVf_3D = np.tile(WfWfVf_3D,(mesh.nz,1,mesh.nx))
        
        UfUf_3D = np.reshape(self.UfUf_spot,(mesh.ny,1))
        UfUf_3D = np.tile(UfUf_3D,(mesh.nz,1,mesh.nx))
        
        UfVf_3D = np.reshape(self.UfVf_spot,(mesh.ny,1))
        UfVf_3D = np.tile(UfVf_3D,(mesh.nz,1,mesh.nx))
        
        VfVf_3D = np.reshape(self.VfVf_spot,(mesh.ny,1))
        VfVf_3D = np.tile(VfVf_3D,(mesh.nz,1,mesh.nx))
        
        WfWf_3D = np.reshape(self.WfWf_spot,(mesh.ny,1))
        WfWf_3D = np.tile(WfWf_3D,(mesh.nz,1,mesh.nx))
        
        
        dUmeandy_3D = D1_3Dy(U_mean_3D * self.spot, mesh.Yc)
        dUmeandy_3D = D1_3Dy_IBM(dUmeandy_3D, U_mean_3D * self.spot, settings, mesh.Yc)

        dUfUfVf_Meandy = D1_3Dy(self.spot * UfUfVf_3D, mesh.Yc)
        dUfUfVf_Meandy = D1_3Dy_IBM(dUfUfVf_Meandy, UfUfVf_3D * self.spot, settings, mesh.Yc)
        
        dUfVfVf_Meandy = D1_3Dy(self.spot * UfVfVf_3D, mesh.Yc)
        dUfVfVf_Meandy = D1_3Dy_IBM(dUfVfVf_Meandy, UfVfVf_3D * self.spot, settings, mesh.Yc)
        
        dVfVfVf_Meandy = D1_3Dy(self.spot * VfVfVf_3D, mesh.Yc)
        dVfVfVf_Meandy = D1_3Dy_IBM(dVfVfVf_Meandy, VfVfVf_3D * self.spot, settings, mesh.Yc)
        
        dWfWfVf_Meandy = D1_3Dy(self.spot * WfWfVf_3D, mesh.Yc)
        dWfWfVf_Meandy = D1_3Dy_IBM(dWfWfVf_Meandy, WfWfVf_3D * self.spot, settings, mesh.Yc)
        
        dUmeandy = average_xz_spot(dUmeandy_3D, self.spot)
        dUfUfVf_Meandy = average_xz_spot(dUfUfVf_Meandy, self.spot)
        dUfVfVf_Meandy = average_xz_spot(dUfVfVf_Meandy, self.spot)
        dVfVfVf_Meandy = average_xz_spot(dVfVfVf_Meandy, self.spot)
        dWfWfVf_Meandy = average_xz_spot(dWfWfVf_Meandy, self.spot)

        d2UfUf_Meandy2 = D2_3Dy(self.spot * UfUf_3D, mesh.Yc, mesh.Y)
        d2UfUf_Meandy2 = D2_3Dy_IBM(d2UfUf_Meandy2, UfUf_3D * self.spot, settings, mesh.Yc, 1/mesh.ny)
        d2UfUf_Meandy2 = average_xz_spot(d2UfUf_Meandy2, self.spot)
        
        d2UfVf_Meandy2 = D2_3Dy(self.spot * UfVf_3D, mesh.Yc, mesh.Y)
        d2UfVf_Meandy2 = D2_3Dy_IBM(d2UfVf_Meandy2, UfVf_3D * self.spot, settings, mesh.Yc, 1/mesh.ny)
        d2UfVf_Meandy2 = average_xz_spot(d2UfVf_Meandy2, self.spot)
        
        d2VfVf_Meandy2 = D2_3Dy(self.spot * VfVf_3D, mesh.Yc, mesh.Y)
        d2VfVf_Meandy2 = D2_3Dy_IBM(d2VfVf_Meandy2, VfVf_3D * self.spot, settings, mesh.Yc, 1/mesh.ny)
        d2VfVf_Meandy2 = average_xz_spot(d2VfVf_Meandy2, self.spot)
        
        d2WfWf_Meandy2 = D2_3Dy(self.spot * WfWf_3D, mesh.Yc, mesh.Y)
        d2WfWf_Meandy2 = D2_3Dy_IBM(d2WfWf_Meandy2, WfWf_3D * self.spot, settings, mesh.Yc, 1/mesh.ny)
        d2WfWf_Meandy2 = average_xz_spot(d2WfWf_Meandy2, self.spot)
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TransportEquationsTerms/'
        createPostProcessingDirectory(current_path_postProcess)

        #dUfUf/dt
        self.uu_budget['Production'] = - 2*self.UfVf_spot * dUmeandy
        self.uu_budget['Turbulent Diffusive Flux'] = - dUfUfVf_Meandy
        self.uu_budget['Molecular diffusion'] = 1/settings.Re * d2UfUf_Meandy2
                
        write_stats(self.uu_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dUfVf/dt
        self.uv_budget['Production'] = - self.VfVf_spot * dUmeandy
        self.uv_budget['Turbulent Diffusive Flux'] = - dUfVfVf_Meandy
        self.uv_budget['Molecular diffusion'] = 1/settings.Re * d2UfVf_Meandy2
                
        write_stats(self.uv_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dVfVf/dt
        self.vv_budget['Turbulent Diffusive Flux'] = - dVfVfVf_Meandy
        self.vv_budget['Molecular diffusion'] = 1/settings.Re * d2VfVf_Meandy2
                
        write_stats(self.vv_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dWfWf/dt
        self.ww_budget['Turbulent Diffusive Flux'] = - dWfWfVf_Meandy
        self.ww_budget['Molecular diffusion'] = 1/settings.Re * d2WfWf_Meandy2
                
        write_stats(self.ww_budget, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def save_transport_equation_terms_ibm(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save velocities budgets (uu,uv,vv,ww) to file.
        Compute reamining quantities needed for the budgets.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        ####################################################
        #Computation of all the gradients needed

        #gradient mean_velocity
        dUmeandy = D1_3Dy(self.U.mean, mesh.Yc)

        dUfUfVf_Meandy = D1_3Dy(self.UfUfVf_spot, mesh.Yc)
        dUfVfVf_Meandy = D1_3Dy(self.UfVfVf_spot, mesh.Yc)
        dVfVfVf_Meandy = D1_3Dy(self.VfVfVf_spot, mesh.Yc)
        dWfWfVf_Meandy = D1_3Dy(self.WfWfVf_spot, mesh.Yc)

        d2UfUf_Meandy2 = D2_3Dy(self.UfUf_spot, mesh.Yc, mesh.Y)
        d2UfVf_Meandy2 = D2_3Dy(self.UfVf_spot, mesh.Yc, mesh.Y)
        d2VfVf_Meandy2 = D2_3Dy(self.VfVf_spot, mesh.Yc, mesh.Y)
        d2WfWf_Meandy2 = D2_3Dy(self.WfWf_spot, mesh.Yc, mesh.Y)
        
        dUmeandy = average_xz_spot(dUmeandy, self.spot)

        dUfUfVf_Meandy = average_xz_spot(dUfUfVf_Meandy, self.spot)
        dUfVfVf_Meandy = average_xz_spot(dUfVfVf_Meandy, self.spot)
        dVfVfVf_Meandy = average_xz_spot(dVfVfVf_Meandy, self.spot)
        dWfWfVf_Meandy = average_xz_spot(dWfWfVf_Meandy, self.spot)

        d2UfUf_Meandy2 = average_xz_spot(d2UfUf_Meandy2, self.spot)
        d2UfVf_Meandy2 = average_xz_spot(d2UfVf_Meandy2, self.spot)
        d2VfVf_Meandy2 = average_xz_spot(d2VfVf_Meandy2, self.spot)
        d2WfWf_Meandy2 = average_xz_spot(d2WfWf_Meandy2, self.spot)
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TransportEquationsTerms/'
        createPostProcessingDirectory(current_path_postProcess)

        #dUfUf/dt
        self.uu_budget['Production'] = - 2*average_xz_spot(self.UfVf_spot, self.spot) * dUmeandy
        self.uu_budget['Turbulent Diffusive Flux'] = - dUfUfVf_Meandy
        self.uu_budget['Molecular diffusion'] = 1/settings.Re * d2UfUf_Meandy2
                
        write_stats(self.uu_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dUfVf/dt
        self.uv_budget['Production'] = - average_xz_spot(self.VfVf_spot, self.spot) * dUmeandy
        self.uv_budget['Turbulent Diffusive Flux'] = - dUfVfVf_Meandy
        self.uv_budget['Molecular diffusion'] = 1/settings.Re * d2UfVf_Meandy2
                
        write_stats(self.uv_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dVfVf/dt
        self.vv_budget['Turbulent Diffusive Flux'] = - dVfVfVf_Meandy
        self.vv_budget['Molecular diffusion'] = 1/settings.Re * d2VfVf_Meandy2
                
        write_stats(self.vv_budget, current_path_postProcess, current_iteration, mesh.ny)

        #dWfWf/dt
        self.ww_budget['Turbulent Diffusive Flux'] = - dWfWfVf_Meandy
        self.ww_budget['Molecular diffusion'] = 1/settings.Re * d2WfWf_Meandy2
                
        write_stats(self.ww_budget, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def save_premultiplied_spectrums(self, settings:Settings, mesh:CFD_mesh, current_iteration = None, TC=False):
        """Save premultiplied spectrums (uu,uv,vv,ww) to file.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        if (TC):
                
            freq_x_theo = (np.arange(self.width_x)/self.Lx_TC)
            freq_z_theo = (np.arange(self.width_z)/self.Lz_TC)
    
            Lx = np.arange(1,np.floor(self.width_x/2),dtype="int")
            Lz = np.arange(1,np.floor(self.width_z/2),dtype="int")
            
            lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
            lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
            
            # store result as individual HDF5 file
            current_path_postProcess = settings.path_postProcess+'/Spectrums/'
            createPostProcessingDirectory(current_path_postProcess)
        
            # read statistic profiles to ascii file
            if (current_iteration is None):
                fnam = 'Spectrums'
            else:
                fnam = 'Spectrums_' + str(current_iteration*1000)
    
            filename = current_path_postProcess + fnam
            
            fields_to_write = {
                'uu_spectrum': self.uu_spectrum,
                'vv_spectrum': self.vv_spectrum,
                'ww_spectrum': self.ww_spectrum,
                'uv_spectrum': self.uv_spectrum,
                'lambda_x_plus': lambda_x_theo[Lx],
                'lambda_z_plus': lambda_z_theo[Lz]
                }
            
            write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)
            
            
        else:
                
            freq_x_theo = (np.arange(mesh.nx)/settings.Lx)
            freq_z_theo = (np.arange(mesh.nz)/settings.Lz)
    
            Lx = mesh.length_x
            Lz = np.arange(1,np.floor(mesh.nz/2),dtype="int")
            
            lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
            lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
            
            # store result as individual HDF5 file
            current_path_postProcess = settings.path_postProcess+'/Spectrums/'
            createPostProcessingDirectory(current_path_postProcess)
        
            # read statistic profiles to ascii file
            if (current_iteration is None):
                fnam = 'Spectrums'
            else:
                fnam = 'Spectrums_' + str(current_iteration*1000)
    
            filename = current_path_postProcess + fnam
            
            fields_to_write = {
                'uu_spectrum': self.uu_spectrum,
                'vv_spectrum': self.vv_spectrum,
                'ww_spectrum': self.ww_spectrum,
                'uv_spectrum': self.uv_spectrum,
                'lambda_x_plus': lambda_x_theo[Lx],
                'lambda_z_plus': lambda_z_theo[Lz]
                }
            
            write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)


    def save_premultiplied_spectrums_ibm(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save premultiplied spectrums (uu,uv,vv,ww) to file while the immersed boundary method is used in the computation.
        Two different files are created, one for 0<y*<k* and one for for k*<y*<1
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """

        # First, write the spectrums for 0<y*<k*
        freq_x_theo = (np.arange(self.width_x)/self.Lx_pattern)
        freq_z_theo = (np.arange(self.width_z)/self.Lz_pattern)

        Lx = np.arange(1,np.floor(self.width_x/2),dtype="int")
        Lz = np.arange(1,np.floor(self.width_z/2),dtype="int")
        
        lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
        lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
        
        # store result as individual HDF5 file
        current_path_postProcess = settings.path_postProcess+'/Spectrums/'
        createPostProcessingDirectory(current_path_postProcess)
    
        # read statistic profiles to ascii file
        if (current_iteration is None):
            fnam = 'Spectrums_pattern'
        else:
            fnam = 'Spectrums_pattern_' + str(current_iteration*1000)

        filename = current_path_postProcess + fnam
        
        fields_to_write = {
            'uu_spectrum': self.uu_spectrum_pattern,
            'vv_spectrum': self.vv_spectrum_pattern,
            'ww_spectrum': self.ww_spectrum_pattern,
            'uv_spectrum': self.uv_spectrum_pattern,
            'lambda_x_plus': lambda_x_theo[Lx],
            'lambda_z_plus': lambda_z_theo[Lz]
            }
        
        write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)
            
            
        # First, write the spectrums for k*<y*<1 
        freq_x_theo = (np.arange(mesh.nx)/settings.Lx)
        freq_z_theo = (np.arange(mesh.nz)/settings.Lz)

        Lx = mesh.length_x
        Lz = np.arange(1,np.floor(mesh.nz/2),dtype="int")
        
        lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
        lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
        
        # store result as individual HDF5 file
        current_path_postProcess = settings.path_postProcess+'/Spectrums/'
        createPostProcessingDirectory(current_path_postProcess)
    
        # read statistic profiles to ascii file
        if (current_iteration is None):
            fnam = 'Spectrums'
        else:
            fnam = 'Spectrums_' + str(current_iteration*1000)

        filename = current_path_postProcess + fnam
        
        fields_to_write = {
            'uu_spectrum': self.uu_spectrum,
            'vv_spectrum': self.vv_spectrum,
            'ww_spectrum': self.ww_spectrum,
            'uv_spectrum': self.uv_spectrum,
            'lambda_x_plus': lambda_x_theo[Lx],
            'lambda_z_plus': lambda_z_theo[Lz]
            }
        
        write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)
        
        
         
    def compute_transport_equation_terms_temperature(self, settings:Settings, mesh:CFD_mesh):
        """Compute temperature budget (1D numpy array with shape (ny)) of the flow in the spot.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        # Compute stats
        
        # dTf/dy in spot
        tmp_dTfdy = D1_3Dy(self.T.get_fluctuations_spot(), mesh.Yc)
        tmp_dTfdz = np.gradient(self.T.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dTfdx = np.gradient(self.T.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
                            
        # dT'²/dx² in the spot
        dTfdx_square_spot = compute_correlations_spot(self.thermal_spot * tmp_dTfdx, tmp_dTfdx)
        
        # dT'²/dy² in the spot
        dTfdy_square_spot = compute_correlations_spot(self.thermal_spot * tmp_dTfdy, tmp_dTfdy)
        
        # dT'²/dz² in the spot
        dTfdz_square_spot = compute_correlations_spot(self.thermal_spot * tmp_dTfdz, tmp_dTfdz)

        #dTfTf/dt
        self.tt_budget['Dissipation'] += - ( 1 / (settings.Pr*settings.Re) ) * (dTfdx_square_spot + dTfdy_square_spot + dTfdz_square_spot)
        
        
         
    def compute_transport_equation_terms_temperature_sedat(self, settings:Settings, mesh:CFD_mesh):
        """Compute temperature budget (1D numpy array with shape (ny)) of the flow in the spot.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
                
        # dTf/dy in spot
        tmp_dTfdy = D1_3Dy(self.thermal_spot * self.T.get_fluctuations_spot(), mesh.Yc)
        tmp_dTfdy = D1_3Dy_IBM(tmp_dTfdy, self.thermal_spot * self.T.get_fluctuations_spot(), settings, mesh.Yc)
        tmp_dTfdz = np.gradient(self.thermal_spot * self.T.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dTfdx = np.gradient(self.thermal_spot * self.T.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
                            
        # dT'²/dx² in the spot
        dTfdx_square_spot = compute_correlations_spot(self.thermal_spot * tmp_dTfdx, tmp_dTfdx)
        
        # dT'²/dy² in the spot
        dTfdy_square_spot = compute_correlations_spot(self.thermal_spot * tmp_dTfdy, tmp_dTfdy)
        
        # dT'²/dz² in the spot
        dTfdz_square_spot = compute_correlations_spot(self.thermal_spot * tmp_dTfdz, tmp_dTfdz)

        #dTfTf/dt
        self.tt_budget['Dissipation'] += - ( 1 / (settings.Pr*settings.Re) ) * (dTfdx_square_spot + dTfdy_square_spot + dTfdz_square_spot)
        
        
    def compute_premultiplied_spectrums_temperature(self, settings:Settings, mesh:CFD_mesh, TC=False):
        """Compute premultiplied spectrums (tt) (3D numpy array with shape (nz,ny/2,nx)) of the flow in the spot.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        if (TC):
                
            # Only use in TC of the spot
            filtre_de_hann_z = np.zeros((mesh.nz,mesh.nx))
            filtre_de_hann_x = np.zeros((mesh.nz,mesh.nx))
            
            shift_index = settings.shift

            spot_centerline = np.roll(self.spot[:,0,:],-shift_index,axis=-1)        

            width_x = np.sum(spot_centerline,axis=1)
            width_z = np.sum(spot_centerline,axis=0)

            max_width_x = np.max(width_x)
            max_width_z = np.max(width_z)
            
            i_start = np.where(spot_centerline[mesh.nz//2,:]==1)[0][0]
            i_end = np.where(spot_centerline[mesh.nz//2,:]==1)[0][-1]
            
            
            k_start = np.where(spot_centerline[:,i_start]==1)[0][0]
            k_end = np.where(spot_centerline[:,i_start]==1)[0][-1]

            for i in range(i_start,i_end+1):
                
                filtre_de_hann_x[:,i] = 0.5 - 0.5 * np.cos(2*np.pi*(i-i_start)/(max_width_x))

            for k in range(k_start,k_end+1):
                
                filtre_de_hann_z[k,:] = 0.5 - 0.5 * np.cos(2*np.pi*(k-k_start)/max_width_z)
                
            window_function = filtre_de_hann_x * filtre_de_hann_z * spot_centerline
            window_function = np.roll(window_function,shift_index,axis=-1)
            
            self.window_function = window_function
            
            # mesh.length_x = np.arange(1,np.floor(max_width_x/2),dtype="int")
            mesh.width_x = max_width_x
        
            # Then compute spectrums ...
            
            # The size of the array for the spectrums has to be redefined
            self.tt_spectrum = np.zeros((int(max_width_z)//2-1,mesh.ny//2,int(max_width_x)//2-1))
            
    
            for j in range(mesh.ny//2):
            
                # T fluctuations
                tmp_slice = np.roll(self.T.get_fluctuations_spot()[:,j,:]*self.window_function,-shift_index,axis=-1)[k_start:k_end+1,i_start:i_end+1]
                fhat_2D_t = np.fft.fft2(tmp_slice)[:,self.length_x]
                fhat_2D_t = fhat_2D_t / (int(max_width_x)*int(max_width_z))
            
                for i in range(int(max_width_x)//2-1):
                    for k in range(int(max_width_z)//2-1):
                        
                        # premultiplied uu spectrum
                        tmp1 = fhat_2D_t[k,i]
                        tmp2 = fhat_2D_t[int(max_width_z)-k-1,i]
                        
                        self.tt_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_TC) * (2*np.pi*i/self.Lx_TC)
    
            del fhat_2D_t
            
            
        else:
        
            for j in range(mesh.ny//2):
            
                # T fluctuations
                fhat_2D_t = np.fft.fft2((self.T.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_t = fhat_2D_t / (mesh.nx*mesh.nz)
            
                for i in range(mesh.nx//2-1):
                    for k in range(mesh.nz//2-1):
                        
                        # premultiplied uu spectrum
                        tmp1 = fhat_2D_t[k,i]
                        tmp2 = fhat_2D_t[mesh.nz-k-1,i]
                        
                        self.tt_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)
            
            del fhat_2D_t
                        
                        
    def compute_premultiplied_spectrums_temperature_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute premultiplied spectrums (tt) (3D numpy array with shape (nz,ny/2,nx)) of the flow while using the immersed boundary method in the domain without the object.
        The idea is apply Hann windowing function for 0<y*<k* on half-pattern and classic FFT computation for k*<y*<1
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        ### Finally, compute the spectrums
        for j in range(mesh.ny//2):

            if (j<self.width_y): # for 0<y*<k*

                for n in range(self.number_of_patterns):

                    # First half of the channel in the spanwise direction
                    i_start_tmp = settings.i_end[int(4*n)] + 1
                    i_end_tmp = i_start_tmp + self.width_x

                    k_start_tmp = settings.k_start[int(4*n)]
                    k_end_tmp = k_start_tmp + self.width_z

                    # u fluctuations
                    fhat_2D_t = np.fft.fft2(self.T.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_t = fhat_2D_t / (int(self.width_x)*int(self.width_z))
                
                    for i in range(int(self.width_x)//2-1):
                        for k in range(int(self.width_z)//2-1):
                            
                            # premultiplied tt spectrum
                            tmp1 = fhat_2D_t[k,i]
                            tmp2 = fhat_2D_t[int(self.width_z)-k-1,i]
                            
                            self.tt_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)


                    # Second half of the channel in the spanwise direction
                    i_start_tmp = settings.i_end[int(4*n+1)] + 1
                    i_end_tmp = i_start_tmp + self.width_x

                    k_start_tmp = settings.k_start[int(4*n+1)]
                    k_end_tmp = k_start_tmp + self.width_z

                    # u fluctuations
                    fhat_2D_t = np.fft.fft2(self.T.get_fluctuations_spot()[k_start_tmp:k_end_tmp,j,i_start_tmp:i_end_tmp]*self.window_function)[:,self.length_x_pattern]
                    fhat_2D_t = fhat_2D_t / (int(self.width_x)*int(self.width_z))
                
                    for i in range(int(self.width_x)//2-1):
                        for k in range(int(self.width_z)//2-1):
                            
                            # premultiplied uu spectrum
                            tmp1 = fhat_2D_t[k,i]
                            tmp2 = fhat_2D_t[int(self.width_z)-k-1,i]
                            
                            self.tt_spectrum_pattern[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/self.Lz_pattern) * (2*np.pi*i/self.Lx_pattern)

            else: # for k*<y*<1

                # u fluctuations
                fhat_2D_t = np.fft.fft2((self.T.get_fluctuations())[:,j,:])[:,mesh.length_x]
                fhat_2D_t = fhat_2D_t / (mesh.nx*mesh.nz)
            
                for i in range(mesh.nx//2-1):
                    for k in range(mesh.nz//2-1):
                        
                        # premultiplied uu spectrum
                        tmp1 = fhat_2D_t[k,i]
                        tmp2 = fhat_2D_t[mesh.nz-k-1,i]
                        
                        self.tt_spectrum[k,j,i] += ( tmp1 * np.conj(tmp1) + tmp2 * np.conj(tmp2) )* (2*np.pi*k/settings.Lz) * (2*np.pi*i/settings.Lx)


        del fhat_2D_t
        
             
    def save_transport_equation_terms_temperature(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save temperature budgets to file.
        Compute reamining quantities needed for the budgets.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # Compute stats
        ####################################################
        tmp_dT_meandy = D1_3Dy(self.T.mean_spot, mesh.Yc)
                                
        # dTf²/dy in spot
        dVfTfTf_spotdy = D1_3Dy(self.VfTfTf_spot, mesh.Yc)
                        
        # d²Tf²/dy² in spot
        d2TfTf_spotdy2 = D2_3Dy(self.TfTf_spot, mesh.Yc, mesh.Y)
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
        createPostProcessingDirectory(current_path_postProcess)

        #dTfTf/dt
        self.tt_budget['Turbulent Transport'] = - 0.5 * dVfTfTf_spotdy
        self.tt_budget['Production'] = - self.VfTf_spot * tmp_dT_meandy
        self.tt_budget['Molecular diffusion'] = ( 1 / (2*settings.Pr*settings.Re) ) * d2TfTf_spotdy2
                        
        write_stats(self.tt_budget, current_path_postProcess, current_iteration, mesh.ny)
        
             
    def save_transport_equation_terms_temperature_sedat(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save temperature budgets to file.
        Compute reamining quantities needed for the budgets.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # Compute stats
        ####################################################
        T_mean_3D = np.reshape(self.T.mean_spot,(mesh.ny,1))
        T_mean_3D = np.tile(T_mean_3D,(mesh.nz,1,mesh.nx))
        
        # remove values in objects
        T_mean_3D = T_mean_3D * self.thermal_spot
        # add constant temperature in objects
        T_mean_3D += -settings.deltaT * (1-self.thermal_spot)
        
        dTmeandy_3D = D1_3Dy(T_mean_3D * self.thermal_spot, mesh.Yc, f_wall=-settings.deltaT)
        dTmeandy_3D = D1_3Dy_IBM(dTmeandy_3D, T_mean_3D * self.thermal_spot, settings, mesh.Yc, f_wall=-settings.deltaT)
        
        dTmeandy = average_xz_spot(dTmeandy_3D, self.thermal_spot)
                                
        # dTf²/dy in spot
        VfTfTf_3D = np.reshape(self.VfTfTf_spot,(mesh.ny,1))
        VfTfTf_3D = np.tile(VfTfTf_3D,(mesh.nz,1,mesh.nx))
        
        dVfTfTf_spotdy_3D = D1_3Dy(VfTfTf_3D * self.thermal_spot, mesh.Yc)
        dVfTfTf_spotdy_3D = D1_3Dy_IBM(dVfTfTf_spotdy_3D, VfTfTf_3D * self.thermal_spot, settings, mesh.Yc)
        
        dVfTfTf_spotdy = average_xz_spot(dVfTfTf_spotdy_3D, self.thermal_spot)
                        
        # d²Tf²/dy² in spot
        TfTf_3D = np.reshape(self.TfTf_spot,(mesh.ny,1))
        TfTf_3D = np.tile(TfTf_3D,(mesh.nz,1,mesh.nx))
        
        d2TfTf_spotdy2 = D2_3Dy(TfTf_3D * self.thermal_spot, mesh.Yc, mesh.Y)
        d2TfTf_spotdy2 = D2_3Dy_IBM(d2TfTf_spotdy2, TfTf_3D * self.thermal_spot, settings, mesh.Yc, 1/mesh.ny)
        d2TfTf_spotdy2 = average_xz_spot(d2TfTf_spotdy2, self.thermal_spot)
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
        createPostProcessingDirectory(current_path_postProcess)

        #dTfTf/dt
        self.tt_budget['Turbulent Transport'] = - 0.5 * dVfTfTf_spotdy
        self.tt_budget['Production'] = - self.VfTf_spot * dTmeandy
        self.tt_budget['Molecular diffusion'] = ( 1 / (2*settings.Pr*settings.Re) ) * d2TfTf_spotdy2
                        
        write_stats(self.tt_budget, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def save_premultiplied_spectrums_temperature(self, settings:Settings, mesh:CFD_mesh, current_iteration = None, TC=False):
        """Save premultiplied spectrums (uu,uv,vv,ww) to file.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        if (TC):
                
            freq_x_theo = (np.arange(self.width_x)/self.Lx_TC)
            freq_z_theo = (np.arange(self.width_z)/self.Lz_TC)
    
            Lx = np.arange(1,np.floor(self.width_x/2),dtype="int")
            Lz = np.arange(1,np.floor(self.width_z/2),dtype="int")
            
            lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
            lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
            
            # store result as individual HDF5 file
            current_path_postProcess = settings.path_postProcess+'/TemperatureSpectrums/'
            createPostProcessingDirectory(current_path_postProcess)
        
            # read statistic profiles to ascii file
            if (current_iteration is None):
                fnam = 'Spectrums'
            else:
                fnam = 'Spectrums_' + str(current_iteration*1000)
    
            filename = current_path_postProcess + fnam
            
            fields_to_write = {
                'tt_spectrum': self.tt_spectrum,
                'lambda_x_plus': lambda_x_theo[Lx],
                'lambda_z_plus': lambda_z_theo[Lz]
                }
            
            write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)
            
            
        else:
                
            freq_x_theo = (np.arange(mesh.nx)/settings.Lx)
            freq_z_theo = (np.arange(mesh.nz)/settings.Lz)
    
            Lx = mesh.length_x
            Lz = np.arange(1,np.floor(mesh.nz/2),dtype="int")
            
            lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
            lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
            
            # store result as individual HDF5 file
            current_path_postProcess = settings.path_postProcess+'/TemperatureSpectrums/'
            createPostProcessingDirectory(current_path_postProcess)
        
            # read statistic profiles to ascii file
            if (current_iteration is None):
                fnam = 'Spectrums'
            else:
                fnam = 'Spectrums_' + str(current_iteration*1000)
    
            filename = current_path_postProcess + fnam
            
            fields_to_write = {
                'tt_spectrum': self.tt_spectrum,
                'lambda_x_plus': lambda_x_theo[Lx],
                'lambda_z_plus': lambda_z_theo[Lz]
                }
            
            write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)
            
            
    def save_premultiplied_spectrums_temperature_ibm(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save premultiplied spectrums (tt) to file while the immersed boundary method is used in the computation.
        Two different files are created, one for 0<y*<k* and one for for k*<y*<1
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # First, write the spectrums for 0<y*<k*
        freq_x_theo = (np.arange(self.width_x)/self.Lx_pattern)
        freq_z_theo = (np.arange(self.width_z)/self.Lz_pattern)

        Lx = np.arange(1,np.floor(self.width_x/2),dtype="int")
        Lz = np.arange(1,np.floor(self.width_z/2),dtype="int")
        
        lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
        lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
        
        # store result as individual HDF5 file
        current_path_postProcess = settings.path_postProcess+'/TemperatureSpectrums/'
        createPostProcessingDirectory(current_path_postProcess)
    
        # read statistic profiles to ascii file
        if (current_iteration is None):
            fnam = 'Spectrums_pattern'
        else:
            fnam = 'Spectrums_pattern_' + str(current_iteration*1000)

        filename = current_path_postProcess + fnam
        
        fields_to_write = {
            'tt_spectrum': self.tt_spectrum_pattern,
            'lambda_x_plus': lambda_x_theo[Lx],
            'lambda_z_plus': lambda_z_theo[Lz]
            }
        
        write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)
            
            
        # First, write the spectrums for k*<y*<1 
        freq_x_theo = (np.arange(mesh.nx)/settings.Lx)
        freq_z_theo = (np.arange(mesh.nz)/settings.Lz)

        Lx = mesh.length_x
        Lz = np.arange(1,np.floor(mesh.nz/2),dtype="int")
        
        lambda_x_theo = (1 / freq_x_theo) * self.Re_tau_spot
        lambda_z_theo = (1 / freq_z_theo) * self.Re_tau_spot
        
        # store result as individual HDF5 file
        current_path_postProcess = settings.path_postProcess+'/TemperatureSpectrums/'
        createPostProcessingDirectory(current_path_postProcess)
    
        # read statistic profiles to ascii file
        if (current_iteration is None):
            fnam = 'Spectrums'
        else:
            fnam = 'Spectrums_' + str(current_iteration*1000)

        filename = current_path_postProcess + fnam
        
        fields_to_write = {
            'tt_spectrum': self.tt_spectrum,
            'lambda_x_plus': lambda_x_theo[Lx],
            'lambda_z_plus': lambda_z_theo[Lz]
            }
        
        write_hdf5(filename, fields_to_write, mesh, settings, Re_tau=1)
        
        
    def save_wall_UQ_stats(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save basic velocity stats (mean, RMS, skewness, flatness and Re_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/VelocityStats/'
        createPostProcessingDirectory(settings.path_postProcess)
        createPostProcessingDirectory(current_path_postProcess)
        
        wall_UQ = {
            'name': 'Wall_Velocities_Stats',
            
            'Re_tau_spot': np.array([self.Re_tau_spot])
            }
        
        write_stats(wall_UQ, current_path_postProcess, current_iteration, 1)
    
    
    def save_drag(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save drag components (total, friction and pressure drag) when the IBM is used ...
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/VelocityStats/'
        createPostProcessingDirectory(settings.path_postProcess)
        createPostProcessingDirectory(current_path_postProcess)
        
        drag = {
            'name': 'drag',
            
            'total_drag': np.array([self.total_drag]),
            'friction_drag': np.array([self.friction_drag]),
            'pressure_drag': np.array([self.pressure_drag])
            }
        
        write_stats(drag, current_path_postProcess, current_iteration, 1)
        
        
    def save_basic_velocity_stats(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save basic velocity stats (mean, RMS, skewness, flatness and Re_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/VelocityStats/'
        createPostProcessingDirectory(settings.path_postProcess)
        createPostProcessingDirectory(current_path_postProcess)
        
        # global Mean_RMS
        self.Mean_RMS = {
            'name': 'Mean_RMS',
            
            'Wall normal coordinate': mesh.Yc,
            'U_mean': self.U.mean_spot,
            'V_mean': self.V.mean_spot,
            'W_mean': self.W.mean_spot,
            'U_rms': np.sqrt(self.UfUf_spot),
            'V_rms': np.sqrt(self.VfVf_spot),
            'W_rms': np.sqrt(self.WfWf_spot)
            }
        
        write_stats(self.Mean_RMS, current_path_postProcess, current_iteration, mesh.ny)
        
        Skewness_Flatness = {
            'name': 'Skewness_Flatness',
            
            'Wall normal coordinate': mesh.Yc,
            'U_skewness': self.UfUfUf_spot/(np.sqrt(self.UfUf_spot)**3),
            'V_skewness': self.VfVfVf_spot/(np.sqrt(self.VfVf_spot)**3),
            'W_skewness': self.WfWfWf_spot/(np.sqrt(self.WfWf_spot)**3),
            'U_flatness': self.UfUfUfUf_spot/(np.sqrt(self.UfUf_spot)**4),
            'V_flatness': self.VfVfVfVf_spot/(np.sqrt(self.VfVf_spot)**4),
            'W_flatness': self.WfWfWfWf_spot/(np.sqrt(self.WfWf_spot)**4)
            }
        
        write_stats(Skewness_Flatness, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def save_local_stats(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save local stats (Re_tau_x, T_tau_x, Nu_x)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/LocalStats/'
        createPostProcessingDirectory(settings.path_postProcess)
        createPostProcessingDirectory(current_path_postProcess)
        
        # global Mean_RMS
        self.Local_stats = {
            'name': 'Local_stats',
            
            'streamwise coordinate': mesh.Xc,
            'Re_tau_x': self.Re_tau_x,
            'T_tau_x': self.T_tau_x,
            'Theta_bulk_x': self.theta_b_x,
            'Nu_x': self.Nu_x,
            'Ek_x': self.Ek_x,
            'TT_x': self.TT_x,
            'Omega_x': self.Omega_x
            }
        
        write_stats(self.Local_stats, current_path_postProcess, current_iteration, mesh.nx)
        
        
    def save_basic_velocity_stats_ibm(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save basic velocity stats (mean, RMS, skewness, flatness and Re_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/VelocityStats/'
        createPostProcessingDirectory(settings.path_postProcess)
        createPostProcessingDirectory(current_path_postProcess)
        
        # global Mean_RMS
        self.Mean_RMS = {
            'name': 'Mean_RMS',
            
            'Wall normal coordinate': mesh.Yc,
            'U_mean': self.U.mean_spot,
            'V_mean': self.V.mean_spot,
            'W_mean': self.W.mean_spot,
            'U_rms': np.sqrt(average_xz_spot(self.UfUf_spot, self.spot)),
            'V_rms': np.sqrt(average_xz_spot(self.VfVf_spot, self.spot)),
            'W_rms': np.sqrt(average_xz_spot(self.WfWf_spot, self.spot))
            }
        
        write_stats(self.Mean_RMS, current_path_postProcess, current_iteration, mesh.ny)
        
        Skewness_Flatness = {
            'name': 'Skewness_Flatness',
            
            'Wall normal coordinate': mesh.Yc,
            'U_skewness': average_xz_spot(self.UfUfUf_spot, self.spot)/(np.sqrt(average_xz_spot(self.UfUf_spot, self.spot)**3)),
            'V_skewness': average_xz_spot(self.VfVfVf_spot, self.spot)/(np.sqrt(average_xz_spot(self.VfVf_spot, self.spot)**3)),
            'W_skewness': average_xz_spot(self.WfWfWf_spot, self.spot)/(np.sqrt(average_xz_spot(self.WfWf_spot, self.spot)**3)),
            'U_flatness': average_xz_spot(self.UfUfUfUf_spot, self.spot)/(np.sqrt(average_xz_spot(self.UfUf_spot, self.spot))**4),
            'V_flatness': average_xz_spot(self.VfVfVfVf_spot, self.spot)/(np.sqrt(average_xz_spot(self.VfVf_spot, self.spot))**4),
            'W_flatness': average_xz_spot(self.WfWfWfWf_spot, self.spot)/(np.sqrt(average_xz_spot(self.WfWf_spot, self.spot))**4)
            }
        
        write_stats(Skewness_Flatness, current_path_postProcess, current_iteration, mesh.ny)
        
        wall_UQ = {
            'name': 'Wall_Velocities_Stats',
            
            'Re_tau_spot': np.array([self.Re_tau_spot])
            }
        
        write_stats(wall_UQ, current_path_postProcess, current_iteration, 1)
        
        
    def save_wall_TQ_stats(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save basic temperature stats (mean, RMS, skewness, flatness, correlations, turbulent Prandtl, Nu and T_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
                
        ####################################################
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
        createPostProcessingDirectory(current_path_postProcess)
                
        wall_TQ = {
            'name': 'Wall_Temperature_Stats',
            
            'T_tau': np.array([self.T_tau_spot]),
            'Nu': np.array([self.Nu_spot])
            }
        
        write_stats(wall_TQ, current_path_postProcess, current_iteration, 1)
        
        
    def save_basic_temperature_stats(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save basic temperature stats (mean, RMS, skewness, flatness, correlations, turbulent Prandtl, Nu and T_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        ####################################################
        # turbulent Prandtl
        if (settings.IBM_flag):
            
            # dU_meandy in flow field
            U_mean_3D = np.reshape(self.U.mean_spot,(mesh.ny,1))
            U_mean_3D = np.tile(U_mean_3D,(mesh.nz,1,mesh.nx))
            
            dUmeandy_3D = D1_3Dy(U_mean_3D * self.spot, mesh.Yc)
            dUmeandy_3D = D1_3Dy_IBM(dUmeandy_3D, U_mean_3D * self.spot, settings, mesh.Yc)
            dU_meandy_spot = average_xz_spot(dUmeandy_3D, self.spot)
                
            # dT_meandy in flow field
            T_mean_3D = np.reshape(self.T.mean_spot,(mesh.ny,1))
            T_mean_3D = np.tile(T_mean_3D,(mesh.nz,1,mesh.nx))
            
            
            # remove values in objects
            T_mean_3D = T_mean_3D * self.thermal_spot
            # add constant temperature in objects
            T_mean_3D += -settings.deltaT * (1-self.thermal_spot)
            
            dTmeandy_3D = D1_3Dy(T_mean_3D * self.thermal_spot, mesh.Yc, f_wall=-settings.deltaT)
            dTmeandy_3D = D1_3Dy_IBM(dTmeandy_3D, T_mean_3D * self.thermal_spot, settings, mesh.Yc, f_wall=-settings.deltaT)
            
            dT_meandy_spot = average_xz_spot(dTmeandy_3D, self.thermal_spot)
            
        else:
            # dU_meandy in the spot
            dU_meandy_spot = D1_3Dy(self.U.mean_spot, mesh.Yc)
            
            # dT_meandy in the spot
            dT_meandy_spot = D1_3Dy(self.T.mean_spot, mesh.Yc, f_wall=-settings.deltaT)
        
        self.turb_Pr = ( -self.UfVf_spot/dU_meandy_spot ) / ( -self.VfTf_spot/dT_meandy_spot )
        
        ####################################################
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
        createPostProcessingDirectory(current_path_postProcess)
        
        Mean_Correlations = {
            'name': 'Mean_Correlations',
            
            'Wall normal coordinate': mesh.Yc,
            'T_mean': self.T.mean_spot,
            'T_rms': np.sqrt(self.TfTf_spot),
            'UfTf': self.UfTf_spot,
            'VfTf': self.VfTf_spot,
            'WfTf': self.WfTf_spot,
            'turb_Pr': self.turb_Pr
            }
        
        write_stats(Mean_Correlations, current_path_postProcess, current_iteration, mesh.ny)
        
        Skewness_Flatness = {
            'name': 'Skewness_Flatness',
            
            'Wall normal coordinate': mesh.Yc,
            'T_skewness': self.TfTfTf_spot/np.sqrt(self.TfTf_spot)**3,
            'T_flatness': self.TfTfTfTf_spot/np.sqrt(self.TfTf_spot)**4
            }
        
        write_stats(Skewness_Flatness, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def save_turb_pr_triple_decomposition(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save basic temperature stats (mean, RMS, skewness, flatness, correlations, turbulent Prandtl, Nu and T_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        ####################################################
        ####################################################
        ## compute what we need
           
        # dU_meandy in flow field
        U_mean_3D = np.reshape(self.U.mean_spot,(mesh.ny,1))
        U_mean_3D = np.tile(U_mean_3D,(mesh.nz,1,mesh.nx))
        
        dUmeandy_3D = D1_3Dy(U_mean_3D * self.spot, mesh.Yc)
        dUmeandy_3D = D1_3Dy_IBM(dUmeandy_3D, U_mean_3D * self.spot, settings, mesh.Yc)
        dU_meandy_spot = average_xz_spot(dUmeandy_3D, self.spot)
            
        # dT_meandy in flow field
        T_mean_3D = np.reshape(self.T.mean_spot,(mesh.ny,1))
        T_mean_3D = np.tile(T_mean_3D,(mesh.nz,1,mesh.nx))
        
        # remove values in objects
        T_mean_3D = T_mean_3D * self.thermal_spot
        # add constant temperature in objects
        T_mean_3D += -settings.deltaT * (1-self.thermal_spot)
        
        dTmeandy_3D = D1_3Dy(T_mean_3D * self.thermal_spot, mesh.Yc, f_wall=-settings.deltaT)
        dTmeandy_3D = D1_3Dy_IBM(dTmeandy_3D, T_mean_3D * self.thermal_spot, settings, mesh.Yc, f_wall=-settings.deltaT)
        
        dT_meandy_spot = average_xz_spot(dTmeandy_3D, self.thermal_spot)
        
        ####################################################
        self.UV_td = average_xz_spot( self.spot * self.U.mean, self.V.mean)
        self.VT_td = average_xz_spot( self.thermal_spot * self.V.mean, self.T.mean)
        
        ####################################################
        
        self.turb_Pr = ( self.UfVf_spot/dU_meandy_spot ) / ( self.VfTf_spot/dT_meandy_spot )
        self.turb_Pr_td = ( self.UV_td/dU_meandy_spot ) / ( self.VT_td/dT_meandy_spot )
        self.turb_Pr_effective = ( (self.UfVf_spot+self.UV_td)/dU_meandy_spot ) / ( (self.VfTf_spot+self.VT_td)/dT_meandy_spot )
        
        ####################################################
        ####################################################
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
        createPostProcessingDirectory(current_path_postProcess)
        
        turb_Pr = {
            'name': 'turb_Pr',
            
            'Wall normal coordinate': mesh.Yc,
            'turb_Pr': self.turb_Pr,
            'turb_Pr_td': self.turb_Pr_td,
            'turb_Pr_effective': self.turb_Pr_effective
            }
        
        write_stats(turb_Pr, current_path_postProcess, current_iteration, mesh.ny)
        
        
        
    def save_basic_temperature_stats_sedat(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save basic temperature stats (mean, RMS, skewness, flatness, correlations, turbulent Prandtl, Nu and T_tau)
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        ####################################################
        # turbulent Prandtl
        
        # dU_meandy in the spot
        dU_meandy_spot = D1_3Dy(self.U.mean_spot, mesh.Yc)
        
        # dT_meandy in the spot
        dT_meandy_spot = D1_3Dy(self.T.mean_spot, mesh.Yc, correct_at_wall = False)
        
        self.turb_Pr = ( -self.UfVf_spot/dU_meandy_spot ) / ( -self.VfTf_spot/dT_meandy_spot )
        
        ####################################################
        # save stats
        current_path_postProcess = settings.path_postProcess+'/TemperatureStats/'
        createPostProcessingDirectory(current_path_postProcess)
        
        Mean_Correlations = {
            'name': 'Mean_Correlations',
            
            'Wall normal coordinate': mesh.Yc,
            'T_mean': self.T.mean_spot,
            'T_rms': np.sqrt(self.TfTf_spot),
            'UfTf': self.UfTf_spot,
            'VfTf': self.VfTf_spot,
            'WfTf': self.WfTf_spot,
            'turb_Pr': self.turb_Pr
            }
        
        write_stats(Mean_Correlations, current_path_postProcess, current_iteration, mesh.ny)
        
        Skewness_Flatness = {
            'name': 'Skewness_Flatness',
            
            'Wall normal coordinate': mesh.Yc,
            'T_skewness': self.TfTfTf_spot/np.sqrt(self.TfTf_spot)**3,
            'T_flatness': self.TfTfTfTf_spot/np.sqrt(self.TfTf_spot)**4
            }
        
        write_stats(Skewness_Flatness, current_path_postProcess, current_iteration, mesh.ny)
        
        wall_TQ = {
            'name': 'Wall_Temperature_Stats',
            
            'T_tau': np.array([self.T_tau_spot]),
            'Nu': np.array([self.Nu_spot])
            }
        
        write_stats(wall_TQ, current_path_postProcess, current_iteration, 1)
        
        
    def save_reynolds_stresses(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save reynolds stresses stats
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/ReynoldsStresses/'
        createPostProcessingDirectory(current_path_postProcess)
        
        Reynolds_stresses = {
            'name': 'Reynolds_stresses',
            
            'Wall normal coordinate': mesh.Yc,
            'uu': self.UfUf_spot,
            'vv': self.VfVf_spot,
            'ww': self.WfWf_spot,
            'uv': self.UfVf_spot,
            'uw': self.UfWf_spot,
            'vw': self.VfWf_spot
            }
        
        write_stats(Reynolds_stresses, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def save_reynolds_stresses_ibm(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save reynolds stresses stats
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/ReynoldsStresses/'
        createPostProcessingDirectory(current_path_postProcess)
        
        Reynolds_stresses = {
            'name': 'Reynolds_stresses',
            
            'Wall normal coordinate': mesh.Yc,
            'uu': average_xz_spot(self.UfUf_spot, self.spot),
            'vv': average_xz_spot(self.VfVf_spot, self.spot),
            'ww': average_xz_spot(self.WfWf_spot, self.spot),
            'uv': average_xz_spot(self.UfVf_spot, self.spot),
            'uw': average_xz_spot(self.UfWf_spot, self.spot),
            'vw': average_xz_spot(self.VfWf_spot, self.spot)
            }
        
        write_stats(Reynolds_stresses, current_path_postProcess, current_iteration, mesh.ny)
        
        
    def compute_vorticity(self, settings:Settings, mesh:CFD_mesh):
        """Compute vorticity stats
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        # compute stats
        ####################################################
        tmp_dUfdz = np.gradient(self.U.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(self.V.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)

        tmp_dVfdx = np.gradient(self.V.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(self.W.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        
        tmp_dUfdy = D1_3Dy(self.U.get_fluctuations_spot(), mesh.Yc)
        tmp_dWfdy = D1_3Dy(self.W.get_fluctuations_spot(), mesh.Yc)
        
        # Omega_x in the spot
        tmp = self.spot * (tmp_dWfdy-tmp_dVfdz)
        tmp[tmp==0] = np.nan
        omega_x_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_y in the spot
        tmp = self.spot * (tmp_dUfdz-tmp_dWfdx)
        tmp[tmp==0] = np.nan
        omega_y_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_z in the spot
        tmp = self.spot * (tmp_dVfdx-tmp_dUfdy)
        tmp[tmp==0] = np.nan
        omega_z_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_x RMS in the spot
        tmp = self.spot * (tmp_dWfdy-tmp_dVfdz)**2
        tmp[tmp==0] = np.nan
        omega_x_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_y RMS in the spot
        tmp = self.spot * (tmp_dUfdz-tmp_dWfdx)**2
        tmp[tmp==0] = np.nan
        omega_y_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_z RMS in the spot
        tmp = self.spot * (tmp_dVfdx-tmp_dUfdy)**2
        tmp[tmp==0] = np.nan
        omega_z_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        del tmp_dUfdz, tmp_dVfdz, tmp_dVfdx, tmp_dWfdx, tmp_dUfdy, tmp_dWfdy, tmp
        
        self.Vort_fluc['omega_x'] += omega_x_spot
        self.Vort_fluc['omega_y'] += omega_y_spot
        self.Vort_fluc['omega_z'] += omega_z_spot
        self.Vort_fluc['omega_x rms'] += np.sqrt(omega_x_rms_spot)
        self.Vort_fluc['omega_y rms'] += np.sqrt(omega_y_rms_spot)
        self.Vort_fluc['omega_z rms'] += np.sqrt(omega_z_rms_spot)

        del omega_x_spot, omega_y_spot, omega_z_spot, omega_x_rms_spot, omega_y_rms_spot, omega_z_rms_spot
        
        
    def compute_vorticity_sedat(self, settings:Settings, mesh:CFD_mesh):
        """Compute vorticity stats
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        # compute stats
        ####################################################
        tmp_dUfdz = np.gradient(self.spot * self.U.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(self.spot * self.V.get_fluctuations_spot(), mesh.Zc, edge_order=2, axis=0)

        tmp_dVfdx = np.gradient(self.spot * self.V.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(self.spot * self.W.get_fluctuations_spot(), mesh.Xc, edge_order=2, axis=2)
        
        tmp_dUfdy = D1_3Dy(self.spot * self.U.get_fluctuations_spot(), mesh.Yc)
        tmp_dWfdy = D1_3Dy(self.spot * self.W.get_fluctuations_spot(), mesh.Yc)
        
        tmp_dUfdy = D1_3Dy_IBM(tmp_dUfdy, self.U.get_fluctuations_spot(), settings, mesh.Yc)
        tmp_dWfdy = D1_3Dy_IBM(tmp_dWfdy, self.W.get_fluctuations_spot(), settings, mesh.Yc)
        
        # Omega_x in the spot
        tmp = self.spot * (tmp_dWfdy-tmp_dVfdz)
        tmp[tmp==0] = np.nan
        omega_x_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_y in the spot
        tmp = self.spot * (tmp_dUfdz-tmp_dWfdx)
        tmp[tmp==0] = np.nan
        omega_y_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_z in the spot
        tmp = self.spot * (tmp_dVfdx-tmp_dUfdy)
        tmp[tmp==0] = np.nan
        omega_z_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_x RMS in the spot
        tmp = self.spot * (tmp_dWfdy-tmp_dVfdz)**2
        tmp[tmp==0] = np.nan
        omega_x_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_y RMS in the spot
        tmp = self.spot * (tmp_dUfdz-tmp_dWfdx)**2
        tmp[tmp==0] = np.nan
        omega_y_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_z RMS in the spot
        tmp = self.spot * (tmp_dVfdx-tmp_dUfdy)**2
        tmp[tmp==0] = np.nan
        omega_z_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        del tmp_dUfdz, tmp_dVfdz, tmp_dVfdx, tmp_dWfdx, tmp_dUfdy, tmp_dWfdy, tmp
        
        self.Vort_fluc['omega_x'] += omega_x_spot
        self.Vort_fluc['omega_y'] += omega_y_spot
        self.Vort_fluc['omega_z'] += omega_z_spot
        self.Vort_fluc['omega_x rms'] += np.sqrt(omega_x_rms_spot)
        self.Vort_fluc['omega_y rms'] += np.sqrt(omega_y_rms_spot)
        self.Vort_fluc['omega_z rms'] += np.sqrt(omega_z_rms_spot)

        del omega_x_spot, omega_y_spot, omega_z_spot, omega_x_rms_spot, omega_y_rms_spot, omega_z_rms_spot
        
        
    def compute_vorticity_ibm(self, settings:Settings, mesh:CFD_mesh):
        """Compute vorticity stats
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        """
        
        # compute stats
        ####################################################
        tmp_dUfdz = np.gradient(self.U.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)
        tmp_dVfdz = np.gradient(self.V.get_fluctuations(), mesh.Zc, edge_order=2, axis=0)

        tmp_dVfdx = np.gradient(self.V.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)
        tmp_dWfdx = np.gradient(self.W.get_fluctuations(), mesh.Xc, edge_order=2, axis=2)
        
        tmp_dUfdy = D1_3Dy(self.U.get_fluctuations(), mesh.Yc)
        tmp_dWfdy = D1_3Dy(self.W.get_fluctuations(), mesh.Yc)
        
        # Omega_x in the spot
        tmp = self.spot * (tmp_dWfdy-tmp_dVfdz)
        tmp[tmp==0] = np.nan
        omega_x_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_y in the spot
        tmp = self.spot * (tmp_dUfdz-tmp_dWfdx)
        tmp[tmp==0] = np.nan
        omega_y_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_z in the spot
        tmp = self.spot * (tmp_dVfdx-tmp_dUfdy)
        tmp[tmp==0] = np.nan
        omega_z_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_x RMS in the spot
        tmp = self.spot * (tmp_dWfdy-tmp_dVfdz)**2
        tmp[tmp==0] = np.nan
        omega_x_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_y RMS in the spot
        tmp = self.spot * (tmp_dUfdz-tmp_dWfdx)**2
        tmp[tmp==0] = np.nan
        omega_y_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        # Omega_z RMS in the spot
        tmp = self.spot * (tmp_dVfdx-tmp_dUfdy)**2
        tmp[tmp==0] = np.nan
        omega_z_rms_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
        
        del tmp_dUfdz, tmp_dVfdz, tmp_dVfdx, tmp_dWfdx, tmp_dUfdy, tmp_dWfdy, tmp
        
        self.Vort_fluc['omega_x'] += omega_x_spot
        self.Vort_fluc['omega_y'] += omega_y_spot
        self.Vort_fluc['omega_z'] += omega_z_spot
        self.Vort_fluc['omega_x rms'] += np.sqrt(omega_x_rms_spot)
        self.Vort_fluc['omega_y rms'] += np.sqrt(omega_y_rms_spot)
        self.Vort_fluc['omega_z rms'] += np.sqrt(omega_z_rms_spot)

        del omega_x_spot, omega_y_spot, omega_z_spot, omega_x_rms_spot, omega_y_rms_spot, omega_z_rms_spot
        
        
    def save_vorticity(self, settings:Settings, mesh:CFD_mesh, current_iteration = None):
        """Save vorticity stats
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
        """
        
        # save stats
        current_path_postProcess = settings.path_postProcess+'/Vorticity/'
        createPostProcessingDirectory(current_path_postProcess)
        
        write_stats(self.Vort_fluc, current_path_postProcess, current_iteration, mesh.ny)
        

    def plot_2D_fields(self, settings:Settings, mesh:CFD_mesh):
        """Plot any 2D fields in any configuration depending on the user settings.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
            
        """
        
        ###########################################################
        #################### W fluctuations
        if (settings.plot_W_fluctuations):
            plot_2D_contours(self.W.get_fluctuations(), 'W fluctuations', mesh, settings)

        #################### V fluctuations
        if (settings.plot_V_fluctuations):
            plot_2D_contours(self.V.get_fluctuations(), 'V fluctuations', mesh, settings) 

        #################### U fluctuations
        if (settings.plot_U_fluctuations):
            plot_2D_contours(self.U.get_fluctuations(), 'U fluctuations', mesh, settings)

        #################### T fluctuations
        if (settings.plot_T_fluctuations):
            plot_2D_contours(self.T.get_fluctuations(), 'T fluctuations', mesh, settings)
            
        ###########################################################
        #################### W instantaneous
        if (settings.plot_W_instantaneous):
            plot_2D_contours(self.W.instantaneous, 'W instantaneous', mesh, settings)
            
        #################### V instantaneous
        if (settings.plot_V_instantaneous):
            plot_2D_contours(self.V.instantaneous, 'V instantaneous', mesh, settings)
            
        #################### U instantaneous
        if (settings.plot_U_instantaneous):
            plot_2D_contours(self.U.instantaneous, 'U instantaneous', mesh, settings)
            
        #################### T instantaneous
        if (settings.plot_T_instantaneous):
            plot_2D_contours(self.T.instantaneous, 'T instantaneous', mesh, settings)
            
            
    def compute_lambda2(self, mesh:CFD_mesh, settings:Settings, Re_tau = 1, current_iteration = None):
        """Compute and save lambda2 to .h5 and xdmf file
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        mesh : CFD_mesh
            Mesh of the simulation that will be written
        Re_tau : Real
            Reynolds shear stress of the simulation. It is used to save in inner (Re_tau != 1) or outer (Re_tau == 1) units
        current_iteration : Integer or None
            Current iteration that will be added to the name of the file
            
        """
        
        # Vertical derivatives
        dUdy = D1_3Dy(self.U.get_fluctuations_spot(), mesh.Yc)*(settings.Re/Re_tau**2)
        dVdy = D1_3Dy(self.V.get_fluctuations_spot(), mesh.Yc)*(settings.Re/Re_tau**2)
        dWdy = D1_3Dy(self.W.get_fluctuations_spot(), mesh.Yc)*(settings.Re/Re_tau**2)
        
        # Streamwise derivative
        dUdx = np.gradient(self.U.get_fluctuations_spot(), mesh.Xc,axis=2)*(settings.Re/Re_tau**2)
        dVdx = np.gradient(self.V.get_fluctuations_spot(), mesh.Xc,axis=2)*(settings.Re/Re_tau**2)
        dWdx = np.gradient(self.W.get_fluctuations_spot(), mesh.Xc,axis=2)*(settings.Re/Re_tau**2)
        
        # Spanwise derivatives
        dUdz = np.gradient(self.U.get_fluctuations_spot(), mesh.Zc,axis=0)*(settings.Re/Re_tau**2)
        dVdz = np.gradient(self.V.get_fluctuations_spot(), mesh.Zc,axis=0)*(settings.Re/Re_tau**2)
        dWdz = np.gradient(self.W.get_fluctuations_spot(), mesh.Zc,axis=0)*(settings.Re/Re_tau**2)
        
        omgX = dWdy - dVdz
        
        print("Computing lambda2...", end='', flush=True)
        tL = timeit.default_timer()
        
        S11 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S11 = dUdx[:,:,:]
        del dUdx
        
        S22 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S22 = dVdy[:,:,:]
        del dVdy
        
        S33 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S33 = dWdz[:,:,:]
        del dWdz
        
        S12 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        W12 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S12 = 0.5*(dUdy[:,:,:] + dVdx[:,:,:])
        W12 = 0.5*(dUdy[:,:,:] - dVdx[:,:,:])
        del dUdy, dVdx
        
        S13 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        W13 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S13 = 0.5*(dUdz[:,:,:] + dWdx[:,:,:])
        W13 = 0.5*(dUdz[:,:,:] - dWdx[:,:,:])
        del dUdz, dWdx
        
        S23 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        W23 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S23 = 0.5*(dVdz[:,:,:] + dWdy[:,:,:])
        W23 = 0.5*(dVdz[:,:,:] - dWdy[:,:,:])
        del dVdz, dWdy
        
        P11 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P12 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P13 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P22 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P23 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P33 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        
        P11 = S11*S11+S12*S12+S13*S13-W12*W12-W13*W13
        P12 = S12*(S11+S22)+S13*S23-W13*W23
        P13 = S13*(S11+S33)+S12*S23+W12*W23
        P22 = S12*S12+S22*S22+S23*S23-W12*W12-W23*W23
        P23 = S23*(S22+S33)+S12*S13-W12*W13
        P33 = S13*S13+S23*S23+S33*S33-W13*W13-W23*W23
        
        del S11, S12, S13, S22, S33
        del W12, W13, W23
        
        aa = -1    
        bb = P11+P22+P33
        cc = P12*P12+P13*P13+P23*P23-P11*P22-P11*P33-P22*P33
        dd = P11*P22*P33+2*P12*P13*P23-P12*P12*P33-P13*P13*P22-P23*P23*P11
        
        del P11, P12, P13, P22, P23, P33
        
        xx = (3*cc/aa-((bb*bb)/(aa*aa)))/3
        yy = ((2*(bb*bb*bb)/(aa*aa*aa)) - (9*(bb*cc)/(aa*aa)) + (27*dd/aa))/27
        zz = ((yy*yy)/4) + ((xx*xx*xx)/27)
        
        pp = bb/(3*aa)
        del aa, bb, cc, dd
        
        ii = np.sqrt((yy*yy/4) - zz)
        jj = -pow(ii,1/3)
        kk = np.arccos(-yy/(2*ii))
        mm = np.cos(kk/3)
        nn = np.sqrt(3)*np.sin(kk/3)
    
        L1 = 2*jj*mm + pp
        L2 = -jj*(mm+nn) + pp
        L3 = -jj*(mm-nn) + pp
        
        del ii, jj, kk, mm, nn
    
        lambda_2 = np.zeros((mesh.nz,mesh.ny,mesh.nx,3))
        lambda_2[:,:,:,0] = L1
        lambda_2[:,:,:,1] = L2
        lambda_2[:,:,:,2] = L3
        
        del L1, L2, L3
        
        lambda_2.sort(axis=-1)
        
        lambda2 = - lambda_2[:,:,:,1]
        
        # S = np.zeros((3,3,mesh.nz,mesh.ny,mesh.nx))
        # A = np.zeros((3,3,mesh.nz,mesh.ny,mesh.nx))
        
        # S[0,0,:,:,:] = dUdx
        # S[0,1,:,:,:] = 0.5*(dUdy + dVdx)
        # S[0,2,:,:,:] = 0.5*(dUdz + dWdx)
        # S[1,1,:,:,:] = dVdy
        # S[1,2,:,:,:] = 0.5*(dVdz + dWdy)
        # S[2,2,:,:,:] = dWdz
        # S[1,0,:,:,:] = S[0,1,:,:,:]
        # S[2,0,:,:,:] = S[0,2,:,:,:]
        # S[2,1,:,:,:] = S[2,1,:,:,:]
        
        # A[0,0,:,:,:] = 0
        # A[0,1,:,:,:] = 0.5*(dUdy - dVdx)
        # A[0,2,:,:,:] = 0.5*(dUdz - dWdx)
        # A[1,1,:,:,:] = 0
        # A[1,2,:,:,:] = 0.5*(dVdz - dWdy)
        # A[2,2,:,:,:] = 0
        # A[1,0,:,:,:] = - A[0,1,:,:,:]
        # A[2,0,:,:,:] = - A[0,2,:,:,:]
        # A[2,1,:,:,:] = - A[2,1,:,:,:]
        
        # del dUdy, dVdy, dWdy
        # del dUdx, dVdx, dWdx
        # del dUdz, dVdz, dWdz
        
        # s_plus_omega = np.zeros((3,3,mesh.nz,mesh.ny,mesh.nx))
        
        # for i in range(3):
        #     for j in range(3):
        #         for k in range(3):
        #             s_plus_omega[i,j,:,:,:] += S[i,k,:,:,:]*S[k,j,:,:,:] + A[i,k,:,:,:]*A[k,j,:,:,:]
                    
        # del A, S
                    
        # eigen_values = np.linalg.eigvals(np.transpose(s_plus_omega))
        
        # eigen_values.sort(axis=-1)
        # lambda2 = np.transpose(eigen_values[:,:,:,1]).real
        
        # del s_plus_omega, eigen_values
        
        print('{:3.1f}'.format(timeit.default_timer()-tL), 'seconds', end='', flush=True)

        print("Saving lambda2...", end='', flush=True)

        y3D = np.reshape(mesh.Yc[:settings.slice_lambda2],(settings.slice_lambda2,1)) * Re_tau
        y3D = np.tile(y3D,(mesh.nz,1,mesh.nx))
        
        # store result as individual HDF5 file
        current_path_postProcess = settings.path_postProcess+'/Lambda2/'
        createPostProcessingDirectory(current_path_postProcess)

        filename = current_path_postProcess + 'Lambda2_' + str(current_iteration*1000)     
        
        print(lambda2.shape)
        
        fields_to_write = {
            'lambda2': lambda2,
            'omgX': omgX,
            'spot': self.spot,
            'y3D': y3D
            }
        
        write_hdf5(filename, fields_to_write, mesh, settings, Re_tau, shift = True)
        
        print('Written file:', filename + '.h5')

        write_xdmf(filename, fields_to_write, mesh, settings)
        
        print('Written file:', filename + '.xmf')

        print("Done!")
