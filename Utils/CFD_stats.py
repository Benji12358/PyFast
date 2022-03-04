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

from stats_utils import compute_correlations_spot, write_stats

from interpolation_schemes import D0s_DRP5_3Dx as interpolate_to_cc_x
from interpolation_schemes import D0s_DRP5_3Dy as interpolate_to_cc_y
from interpolation_schemes import D0s_DRP5_3Dz as interpolate_to_cc_z

from derivative_schemes import D1_3Dy, D2_3Dy

from utilities import createPostProcessingDirectory, find_nearest
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
        
        mean_folder = settings.root + '/' + settings.current_path + '/' + settings.current_path + '/Results/3D/field' + str(settings.preliminary_iterations*1000)

        # U (streamwise velocity here)
        fu = h5py.File(mean_folder + '/W.h5', 'r')
        U_mean = np.array(fu['W'])[:-1,:-1,:-1]
        fu.close()
        # V
        fv = h5py.File(mean_folder + '/V.h5', 'r')
        V_mean = np.array(fv['V'])[:-1,:-1,:-1]
        fv.close()
        # U (spanwise velocity here)
        fw = h5py.File(mean_folder + '/U.h5', 'r')
        W_mean = np.array(fw['U'])[:-1,:-1,:-1]
        fw.close()
        # Pressure
        fp = h5py.File(mean_folder + '/P.h5', 'r')
        P_mean = np.array(fp['P'])[:-1,:-1,:-1]
        fp.close()
        
        # from staggered to cell center field
        U_mean = interpolate_to_cc_x(U_mean)[:,:,:mesh.nx]
        V_mean = interpolate_to_cc_y(V_mean)[:,:,:mesh.nx]
        W_mean = interpolate_to_cc_z(W_mean)[:,:,:mesh.nx]
        
        if (settings.get_T):
            ft = h5py.File(mean_folder + '/sca1.h5', 'r')
            T_mean = np.array(ft['sca1'])[:mesh.nz,:mesh.ny,:mesh.nx]
            ft.close()
        
        print(mean_folder)
        
        if (settings.current_path in ['PL_Vortices_LargeEps_fine', 'PL_Vortices_LargeEps_Ka']):
            
            U_backup = U_mean[mesh.nz//2-1,:,mesh.nx//2-1]
            P_backup = P_mean[mesh.nz//2-1,:,mesh.nx//2-1]
                                            
            W_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            V_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            U_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            P_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
            for i in range(mesh.nx):
                for k in range(mesh.nz):
                    U_mean[k,:,i] = U_backup
                    P_mean[k,:,i] = P_backup
                    
            if (settings.get_T):
                
                T_backup = T_mean[mesh.nz//2-1,:,mesh.nx//2-1]
                
                T_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
                for i in range(mesh.nx):
                    for k in range(mesh.nz):
                        T_mean[k,:,i] = T_backup
                        
        self.U.set_mean(U_mean)
        self.V.set_mean(V_mean)
        self.W.set_mean(W_mean)
        self.P.set_mean(P_mean[:,:,:mesh.nx])
            
        del U_mean, V_mean, W_mean, P_mean
            
        if ('Ka' in settings.current_path) and (settings.get_T):
                
            T_mean = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        
            for j in range(mesh.ny):
                    T_mean[:,j,:] = settings.Pr * np.sqrt(2*settings.Re) * (- 0.5*mesh.Yc[j]**3 + 1/8*mesh.Yc[j]**4 + mesh.Yc[j])
        
        if (settings.get_T):
                            
            self.T.set_mean(T_mean)
                
            del T_mean
            
            
    def compute_mean_fields(self, settings:Settings, mesh:CFD_mesh):
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
            
            mean_folder = settings.root + '/' + settings.current_path + '/' + settings.current_path + '/Results/3D/field' + str(i*1000)
    
            # U (streamwise velocity here)
            fu = h5py.File(mean_folder + '/W.h5', 'r')
            U_tmp = np.array(fu['W'])[:-1,:-1,:-1]
            fu.close()
            # V
            fv = h5py.File(mean_folder + '/V.h5', 'r')
            V_tmp = np.array(fv['V'])[:-1,:-1,:-1]
            fv.close()
            # U (spanwise velocity here)
            fw = h5py.File(mean_folder + '/U.h5', 'r')
            W_tmp = np.array(fw['U'])[:-1,:-1,:-1]
            fw.close()
            # Pressure
            fp = h5py.File(mean_folder + '/P.h5', 'r')
            P_mean += np.array(fp['P'])[:-1,:-1,:-1]
            fp.close()
            
            # from staggered to cell center field
            U_mean += interpolate_to_cc_x(U_tmp)[:,:,:mesh.nx]
            V_mean += interpolate_to_cc_y(V_tmp)[:,:,:mesh.nx]
            W_mean += interpolate_to_cc_z(W_tmp)[:,:,:mesh.nx]
            
            if (settings.get_T):
                ft = h5py.File(mean_folder + '/sca1.h5', 'r')
                T_mean += np.array(ft['sca1'])[:mesh.nz,:mesh.ny,:mesh.nx]
                ft.close()
        
        self.U.set_mean(U_mean/len(settings.all_iterations))
        self.V.set_mean(V_mean/len(settings.all_iterations))
        self.W.set_mean(W_mean/len(settings.all_iterations))
        self.P.set_mean(P_mean[:,:,:mesh.nx]/len(settings.all_iterations))
            
        del U_mean, V_mean, W_mean, P_mean
        
        # Re_tau in the spot
        tmp = self.spot[:,0,:] * self.U.mean[:,0,:]
        tmp[tmp==0] = np.nan
        self.Re_tau_spot = np.sqrt( settings.Re * np.nanmean(tmp) / mesh.Yc[0])
        
        if (settings.symmetry):
            tmp = self.spot[:,-1,:] * self.U.mean[:,-1,:]
            tmp[tmp==0] = np.nan
            self.Re_tau_spot += np.sqrt( settings.Re * np.nanmean(tmp) / mesh.Yc[0])
            self.Re_tau_spot = self.Re_tau_spot/2

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
        
        del U_mean_spot, V_mean_spot, W_mean_spot, P_mean_spot
        
        if (settings.get_T):
                            
            self.T.set_mean(T_mean/len(settings.all_iterations))
                
            del T_mean
            
            # T_tau in the spot
            tmp = self.thermal_spot[:,0,:] * self.T.mean[:,0,:]
            tmp[tmp==0] = np.nan
            self.T_tau_spot = abs((np.nanmean(tmp)+settings.deltaT) / (mesh.Yc[0]*self.Re_tau_spot*settings.Pr))
        
            if (settings.symmetry):
                tmp = self.thermal_spot[:,-1,:] * self.T.mean[:,-1,:]
                tmp[tmp==0] = np.nan
                self.T_tau_spot += abs((np.nanmean(tmp)-settings.deltaT) / (mesh.Yc[0]*self.Re_tau_spot*settings.Pr))
                self.T_tau_spot = self.T_tau_spot/2
                
            # Nusselt in the spot
            tmp = self.thermal_spot[:,0,:] * self.T.mean[:,0,:]
            tmp[tmp==0] = np.nan
            self.Nu_spot = abs((np.nanmean(tmp)+settings.deltaT) / (mesh.Yc[0])) * ( 1 / settings.deltaT )
        
            if (settings.symmetry):
                tmp = self.thermal_spot[:,-1,:] * self.T.mean[:,-1,:]
                tmp[tmp==0] = np.nan
                self.Nu_spot += abs((np.nanmean(tmp)-settings.deltaT) / (mesh.Yc[0])) * ( 1 / settings.deltaT )
                self.Nu_spot = self.Nu_spot/2
                
            ####################################################
            # T_mean in the spot
            tmp = self.thermal_spot * self.T.mean
            tmp[tmp==0] = np.nan
            T_mean_spot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
    
            self.T.set_mean_spot(T_mean_spot)
            
            del T_mean_spot
        

    def read_current_fields(self, settings:Settings, current_iteration):
        """Read the instantaneous fields corresponding the current_iteration of the simulation.
        
        Parameters
        ----------
        settings : Settings
            Settings of the simulation
        current_iteration : Integer
            Iteration corresponding to the fields read and processed
    
            Set instantaneous fields of the CFD_stats object
        """
        
        results_path = settings.root + '/' + settings.current_path + '/' + settings.current_path + '/Results/3D/field' + str(current_iteration*1000)
        
        progress = round(100 * (np.where(settings.all_iterations==current_iteration)[0][0] + 1) / len(settings.all_iterations), 2)
        print(str(progress) + " %")
        
        # U (streamwise velocity here)
        fu = h5py.File(results_path + '/W.h5', 'r')
        U_instantaneous = np.array(fu['W'])[:-1,:-1,:-1]
        fu.close()
        # V
        fv = h5py.File(results_path + '/V.h5', 'r')
        V_instantaneous = np.array(fv['V'])[:-1,:-1,:-1]
        fv.close()
        # W (spanwise velocity here)
        fw = h5py.File(results_path + '/U.h5', 'r')
        W_instantaneous = np.array(fw['U'])[:-1,:-1,:-1]
        fw.close()
        # Pressure
        fp = h5py.File(results_path + '/P.h5', 'r')
        P_instantaneous = np.array(fp['P'])[:-1,:-1,:-1]
        fp.close()
        
        if (settings.get_T):
            ft = h5py.File(results_path + '/sca1.h5', 'r')
            self.T.set_instantaneous(np.array(ft['sca1'])[:-1,:-1,:-1])
            ft.close()
            
        # from staggered to cell center field
        U_instantaneous = interpolate_to_cc_x(U_instantaneous)
        V_instantaneous = interpolate_to_cc_y(V_instantaneous)
        W_instantaneous = interpolate_to_cc_z(W_instantaneous)
        
        self.U.set_instantaneous(U_instantaneous)
        self.V.set_instantaneous(V_instantaneous)
        self.W.set_instantaneous(W_instantaneous)
        self.P.set_instantaneous(P_instantaneous)
        
        del U_instantaneous, V_instantaneous, W_instantaneous, P_instantaneous
        
        
    def get_turbulent_spot(self, settings:Settings, mesh:CFD_mesh):
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
        E_fluctuations = gaussian_filter(E_fluctuations, 6)
        self.spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))
                    
        for j in range(mesh.ny):
            E_max = np.max(E_fluctuations[:,j,:])
            
            alpha = 0.05
            self.spot[:,j,:] = (E_fluctuations[:,j,:]>=alpha*E_max)
            
        ##################################################################
        # Get a more accurate spot
        ##################################################################

        if (settings.heart_spot or settings.wave_packet_spot):
        
            noisy = self.U.get_fluctuations()[:,0,:]**2 + self.V.get_fluctuations()[:,0,:]**2 + self.W.get_fluctuations()[:,0,:]**2
            
            freq_x = np.fft.fftfreq(mesh.nx,mesh.xstep)
            freq_z = np.fft.fftfreq(mesh.nz,mesh.zstep)
            
            lambda_x = 2*np.pi/(0.7*10**(0)) # filter all scale greater than h
            lambda_z = 2*np.pi/(0.7*10**(0)) # filter all scale greater than h
            
            gx = np.exp(-((2*np.pi*freq_x)**2*lambda_x**2/24))      # Gauss exp kernel in x direction
            gz = np.exp(-((2*np.pi*freq_z )**2*lambda_z**2/24))     # Gauss exp kernel in z direction
            g2d = np.outer(gz, gx)                                  # 2d filter kernel in x-z plane
            
            filtered=np.fft.ifft2(np.fft.fft2(noisy)*g2d)
                        
            wall_spot = np.zeros((mesh.nz,mesh.nx))
                                                                    
            alpha = 0.6
            wall_spot[:,:] = (filtered[:,:]>=alpha*np.max(filtered))

            if (settings.heart_spot):
            
                for j in range(mesh.ny):
                    self.spot[:,j,:] = wall_spot[:,:]*self.spot[:,j,:]

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
            
    
    def get_turbulent_thermal_spot(self, settings:Settings, mesh:CFD_mesh):
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
        
        T_fluctuations_spot = gaussian_filter(self.T.get_fluctuations()**2, 6)

        self.thermal_spot = np.zeros((mesh.nz,mesh.ny,mesh.nx))
            
        alpha = 0.05

        for j in range(mesh.ny):
            T_max = np.max(T_fluctuations_spot[:,j,:])
            self.thermal_spot[:,j,:] = (T_fluctuations_spot[:,j,:]>=alpha*T_max)
            
        ##################################################################
        # Get a more accurate spot
        ##################################################################

        if (settings.heart_spot or settings.wave_packet_spot):
        
            noisy = self.T.get_fluctuations()[:,0,:]**2
            
            freq_x = np.fft.fftfreq(mesh.nx,mesh.xstep)
            freq_z = np.fft.fftfreq(mesh.nz,mesh.zstep)
            
            lambda_x = 2*np.pi/(0.7*10**(0)) # filter all scale greater than h
            lambda_z = 2*np.pi/(0.7*10**(0)) # filter all scale greater than h
            
            gx = np.exp(-((2*np.pi*freq_x)**2*lambda_x**2/24))      # Gauss exp kernel in x direction
            gz = np.exp(-((2*np.pi*freq_z )**2*lambda_z**2/24))     # Gauss exp kernel in z direction
            g2d = np.outer(gz, gx)                                  # 2d filter kernel in x-z plane
            
            filtered=np.fft.ifft2(np.fft.fft2(noisy)*g2d)
                        
            wall_spot = np.zeros((mesh.nz,mesh.nx))
                                                                    
            alpha = 0.6
            wall_spot[:,:] = (filtered[:,:]>=alpha*np.max(filtered))

            if (settings.heart_spot):
            
                for j in range(mesh.ny):
                    self.thermal_spot[:,j,:] = wall_spot[:,:]*self.thermal_spot[:,j,:]

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
        
         
    def compute_transport_equation_terms_temperature(self, settings:Settings, mesh:CFD_mesh:
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
            'U_skewness': self.UfUfUf_spot/(np.sqrt(self.UfUf_spot**3)),
            'V_skewness': self.VfVfVf_spot/(np.sqrt(self.VfVf_spot**3)),
            'W_skewness': self.WfWfWf_spot/(np.sqrt(self.WfWf_spot**3)),
            'U_flatness': self.UfUfUfUf_spot/(np.sqrt(self.UfUf_spot)**4),
            'V_flatness': self.VfVfVfVf_spot/(np.sqrt(self.VfVf_spot)**4),
            'W_flatness': self.WfWfWfWf_spot/(np.sqrt(self.WfWf_spot)**4)
            }
        
        write_stats(Skewness_Flatness, current_path_postProcess, current_iteration, mesh.ny)
        
        wall_UQ = {
            'name': 'Wall_Velocities_Stats',
            
            'Re_tau_spot': np.array([self.Re_tau_spot])
            }
        
        write_stats(wall_UQ, current_path_postProcess, current_iteration, 1)
        
        
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
            
            
    def compute_lambda2(self, settings:Settings, mesh:CFD_mesh, Re_tau = 1, current_iteration = None):
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
        dUdy = D1_3Dy(self.U.instantaneous, mesh.Yc)*(settings.Re/Re_tau**2)
        dVdy = D1_3Dy(self.V.instantaneous, mesh.Yc)*(settings.Re/Re_tau**2)
        dWdy = D1_3Dy(self.W.instantaneous, mesh.Yc)*(settings.Re/Re_tau**2)
        
        # Streamwise derivative
        dUdx = np.gradient(self.U.instantaneous, mesh.Xc,axis=2)*(settings.Re/Re_tau**2)
        dVdx = np.gradient(self.V.instantaneous, mesh.Xc,axis=2)*(settings.Re/Re_tau**2)
        dWdx = np.gradient(self.W.instantaneous, mesh.Xc,axis=2)*(settings.Re/Re_tau**2)
        
        # Spanwise derivatives
        dUdz = np.gradient(self.U.instantaneous, mesh.Zc,axis=0)*(settings.Re/Re_tau**2)
        dVdz = np.gradient(self.V.instantaneous, mesh.Zc,axis=0)*(settings.Re/Re_tau**2)
        dWdz = np.gradient(self.W.instantaneous, mesh.Zc,axis=0)*(settings.Re/Re_tau**2)
        
        print("Computing lambda2...", end='', flush=True)
        tL = timeit.default_timer()  
                
        S11 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S12 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S13 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S22 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S23 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        S33 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        
        W12 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        W13 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        W23 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        
        P11 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P12 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P13 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P22 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P23 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        P33 = np.zeros((mesh.nz,mesh.ny,mesh.nx))
        
        S11 = dUdx
        S12 = 0.5*(dUdy + dVdx)
        S13 = 0.5*(dUdz + dWdx)
        S22 = dVdy
        S23 = 0.5*(dVdz + dWdy)
        S33 = dWdz
        
        W12 = 0.5*(dUdy - dVdx)
        W13 = 0.5*(dUdz - dWdx)
        W23 = 0.5*(dVdz - dWdy)
        
        
        P11 = S11*S11+S12*S12+S13*S13-W12*W12-W13*W13
        P12 = S12*(S11+S22)+S13*S23-W13*W23
        P13 = S13*(S11+S33)+S12*S23+W12*W23
        P22 = S12*S12+S22*S22+S23*S23-W12*W12-W23*W23
        P23 = S23*(S22+S33)+S12*S13-W12*W13
        P33 = S13*S13+S23*S23+S33*S33-W13*W13-W23*W23
        
        s_plus_omega = np.array( [[P11,P12,P13],
                               [P12,P22,P23],
                               [P13,P23,P33]] )
        
        del dUdy, dVdy, dWdy
        del dUdx, dVdx, dWdx
        del dUdz, dVdz, dWdz
        
        del S11, S12, S13, S22, S23, S33
        del W12, W13, W23
        del P11, P12, P13, P22, P23, P33
                    
        eigen_values = np.linalg.eigvals(np.transpose(s_plus_omega))
        
        eigen_values.sort(axis=-1)
        lambda2 = np.transpose(eigen_values[:,:,:,1])
        
        del s_plus_omega, eigen_values
        
        print('{:3.1f}'.format(timeit.default_timer()-tL), 'seconds', end='', flush=True)

        print("Saving lambda2...", end='', flush=True)

        y3D = np.reshape(mesh.Yc[:settings.slice_lambda2],(settings.slice_lambda2,1))
        y3D = np.tile(y3D,(mesh.nz,1,mesh.nx))
        
        # store result as individual HDF5 file
        current_path_postProcess = settings.path_postProcess+'/Lambda2/'
        createPostProcessingDirectory(current_path_postProcess)

        filename = current_path_postProcess + 'Lambda2_' + str(current_iteration*1000)     
        
        print(lambda2.shape)
        
        fields_to_write = {
            'lambda2': lambda2,
            'y3D': y3D
            }
        
        write_hdf5(filename, fields_to_write, mesh, settings, Re_tau)
        
        print('Written file:', filename + '.h5')

        write_xdmf(filename, fields_to_write, mesh, settings)
        
        print('Written file:', filename + '.xmf')

        print("Done!")
