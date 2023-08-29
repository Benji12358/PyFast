"""
Created on Mon Feb 28 09:40:01 2022

@author: Benjamin Arrondeau

@title: Utilities functions for statistics

@description: Contains general functions for computing/writing/reading statistics.
"""

import numpy as np
import h5py
import xml.etree.cElementTree as ET
from xml.dom import minidom
import os

from utilities import find_nearest

from CFD_mesh import CFD_mesh
from Settings import Settings

def compute_correlations_spot(array1, array2, array3 = None, array4 = None):
    """Computation of the correlations in the "spot/domain".
    
    Depending on the inputs, compute the 2nd, 3rd of 4th order correlation.
    By default, computes 2nd order correlation.
    
    Parameters
    ----------
    array1 : 3D numpy array with shape (nz,ny,nx)
        A flow field
    array2 : 3D numpy array with shape (nz,ny,nx)
        A flow field
    array3 : 3D numpy array with shape (nz,ny,nx)
        A flow field. Default value of None
    array4 : 3D numpy array with shape (nz,ny,nx)
        A flow field. Default value of None

    Returns
    -------
    1D numpy array with shape (ny)
        The spatial averaged of the computed correlation.
    """
    
    if (array4 is not None):
        
        tmp = array1 * array2 * array3 * array4
        
    elif (array3 is not None):
        
        tmp = array1 * array2 * array3
        
    else:
        
        tmp = array1 * array2
        
    tmp[tmp==0] = np.nan
    correlSpot = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
    
    return correlSpot

def compute_correlations_x(array1, array2, array3 = None, array4 = None):
    """Computation of the correlations with respect to x.
    
    Depending on the inputs, compute the 2nd, 3rd of 4th order correlation.
    By default, computes 2nd order correlation.
    
    Parameters
    ----------
    array1 : 3D numpy array with shape (nz,ny,nx)
        A flow field
    array2 : 3D numpy array with shape (nz,ny,nx)
        A flow field
    array3 : 3D numpy array with shape (nz,ny,nx)
        A flow field. Default value of None
    array4 : 3D numpy array with shape (nz,ny,nx)
        A flow field. Default value of None

    Returns
    -------
    1D numpy array with shape (ny)
        The spatial averaged of the computed correlation.
    """
    
    if (array4 is not None):
        
        tmp = array1 * array2 * array3 * array4
        
    elif (array3 is not None):
        
        tmp = array1 * array2 * array3
        
    else:
        
        tmp = array1 * array2
        
    tmp[tmp==0] = np.nan
    correlx = np.nan_to_num(np.nanmean(tmp,axis=(0)))
    
    return correlx

def compute_correlations_spot_3D(array1, array2, array3 = None, array4 = None):
    """Computation of the correlations.
    
    Depending on the inputs, compute the 2nd, 3rd of 4th order correlation.
    By default, computes 2nd order correlation.
    
    Parameters
    ----------
    array1 : 3D numpy array with shape (nz,ny,nx)
        A flow field
    array2 : 3D numpy array with shape (nz,ny,nx)
        A flow field
    array3 : 3D numpy array with shape (nz,ny,nx)
        A flow field. Default value of None
    array4 : 3D numpy array with shape (nz,ny,nx)
        A flow field. Default value of None

    Returns
    -------
    3D numpy array with shape (nz,ny,nx)
        The spatial averaged of the computed correlation.
    """
    
    if (array4 is not None):
        
        correlSpot = array1 * array2 * array3 * array4
        
    elif (array3 is not None):
        
        correlSpot = array1 * array2 * array3
        
    else:
        
        correlSpot = array1 * array2
    
    return correlSpot

def average_xz_spot(array1, spot):
    """Computation of the average on planes within a spot.
    
    Parameters
    ----------
    array1 : 3D numpy array with shape (nz,ny,nx)
        A flow field
    spot : 3D numpy array with shape (nz,ny,nx)
        The spot/sub domain where the average is computed

    Returns
    -------
    1D numpy array with shape (ny)
        The spatial averaged of the computed correlation.
    """
    
    
    tmp = array1 * spot
    tmp[tmp==0] = np.nan
    avg_array = np.nan_to_num(np.nanmean(tmp,axis=(0,2)))
    
    return avg_array

def write_stats(dictionnary, current_path, current_iteration, ny):
    """Write statistics to .dat file.
    
    Write a dictionnary in file.
    
    Parameters
    ----------
    dictionnary : Dictionnary
        A dictionnary containing all variables to write in the file
    current_path : String
        Path where to write the file
    current_iteration : Integer
        Current iteration that will be added to the name of the file
    ny : Integer
        Number of points to write in file
        
        The dictionnary is formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}
            
    """
    
    os.chdir(current_path)
    
    # write statistic profiles to ascii file
    if (current_iteration is None):
        file = open(str(dictionnary['name']) + '.dat', "w")
    else:
        file = open(str(dictionnary['name']) + str(current_iteration*1000) + '.dat', "w")
    
    counter = 1
    
    for key, array in dictionnary.items():
        
        if key=="name":
            file.write("# column n° " + str(counter) + ": " + str(array) + "\n")
        else:
            file.write("# column n° " + str(counter) + ": " + str(key) + "\n")
        counter += 1
    
    for jth in range(ny):
        
        for key, array in dictionnary.items():
            
            if key!="name":
                file.write(" %23.16e" % (array[jth]))
            
        file.write("\n")   
        
    file.close()
    print("Written " + str(dictionnary['name']) + " to file")

def read_stats(dictionnary, current_path, current_iteration = None):
    """Read statistics contained in .dat file.
    
    Read a file and store value in a dictionnary
    
    Parameters
    ----------
    dictionnary : Dictionnary
        A dictionnary containing all variables to read in the file
    current_path : String
        Path where to read the file
    current_iteration : Integer
        Current iteration that will be added to the name of the file
        
        The dictionnary is formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'variable_1': array of shape (ny),
                ...
                'variable_n': array of shape (ny)}

    Returns
    -------
    Dictionnary
        The dictionnary with all arrays updated to the values contained in the file
            
    """
    
    os.chdir(current_path)
    
    # read statistic profiles to ascii file
    if (current_iteration is None):
        fnam = str(dictionnary['name']) + '.dat'
    else:
        fnam = str(dictionnary['name']) + str(current_iteration*1000) + '.dat'
        
    counter = 0
    number_rows = len(dictionnary.keys())
    # print(fnam)
    
    for key in dictionnary.keys():
        
        if key!="name":
            dictionnary[key] = np.loadtxt(fnam,usecols=counter,skiprows=number_rows)
            counter += 1
    
    print("Read " + str(dictionnary['name']) + " file")
    
    return dictionnary


def read_U_main_quantities(current_path, current_iteration = None):
    """Read main velocity statistics contained in .dat file.
    
    Read the Re_tau of the simulation
    
    Parameters
    ----------
    current_path : String
        Path where to read the file
    current_iteration : Integer
        Current iteration that will be added to the name of the file
        
        The dictionnary is formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'Re_tau_spot': array of shape (1)}

    Returns
    -------
    Dictionnary
        The dictionnary with main velocity statistics updated to the values contained in the file
            
    """
    
    main_ui_quantity = {
        'name': 'Wall_Velocities_Stats',
        
        'Re_tau_spot': None
        }
    
    main_ui_quantity = read_stats(main_ui_quantity, current_path, current_iteration)
    
    return main_ui_quantity


def read_T_main_quantities(current_path, current_iteration = None):
    """Read main temperature statistics contained in .dat file.
    
    Read the T_tau and Nu of the simulation
    
    Parameters
    ----------
    current_path : String
        Path where to read the file
    current_iteration : Integer
        Current iteration that will be added to the name of the file
        
        The dictionnary is formatted like this:
            
            dictionnary = { 
                'name': 'name_of_the_file',
                'T_tau_spot': array of shape (1),
                'Nu_tau_spot': array of shape (1)}

    Returns
    -------
    Dictionnary
        The dictionnary with main temperature statistics updated to the values contained in the file
            
    """
    
    main_T_quantity = {
        'name': 'Wall_Temperature_Stats',
        
        'T_tau_spot': None,
        'Nu_spot': None
        }
    
    main_T_quantity = read_stats(main_T_quantity, current_path, current_iteration)
    
    return main_T_quantity
    
    
def prettify(elem):
    # """Return a pretty-printed XML string for the Element."""
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
    

def write_hdf5(fname, dictionnary, mesh:CFD_mesh, settings:Settings, Re_tau = 1, shift = False):
    """Write 3D fields to .h5 file.
    
    Write a dictionnary of 3D fields in .h5 file.
    For each file, add the mesh to the .h5 file so that it can be read by Paraview.
    
    Parameters
    ----------
    fname : String
        Name of the file to write
    dictionnary : Dictionnary
        A dictionnary containing all variables to write in the file
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    settings : Settings
        Settings of the simulation
    Re_tau : Integer
        Reynolds shear stress of the simulation. It is used to save in inner (Re_tau != 1) or outer (Re_tau == 1) units
        
        The dictionnary is formatted like this:
            
            dictionnary = {
                'variable_1': array of shape (nz,ny,nx),
                ...
                'variable_n': array of shape (nz,ny,nx)}
            
    """
    
    out = h5py.File(fname + '.h5', 'w') # open HDF5 file for writing

    for key, array in dictionnary.items():
        if (len(array.shape)==3):
            if shift:
                out[str(key)] = np.roll(array[:,:settings.slice_lambda2,:],-settings.shift)
            else:
                out[str(key)] = array[:,:settings.slice_lambda2,:]
        else:
            out[str(key)] = array

    # write mesh in h5 file
    # out['Xc'] = (mesh.Xc + mesh.Xc[find_nearest(mesh.Xc,settings.shift)]) * Re_tau
    out['Xc'] = mesh.Xc * Re_tau
    out['Yc'] = mesh.Yc[:settings.slice_lambda2] * Re_tau
    out['Zc'] = mesh.Zc * Re_tau

    out.close() # close HDF5 file
    

def write_xdmf(fname, dictionnary, mesh:CFD_mesh, settings:Settings):
    """Write xdmf file.
    
    Write the corresponded xdmf file to the .h5 file written
    For each file, declare the mesh so that it can be read by Paraview.
    
    Parameters
    ----------
    fname : String
        Name of the file to write
    dictionnary : Dictionnary
        A dictionnary containing all variables to write in the file
    mesh : CFD_mesh
        Mesh of the simulation that will be written
    settings : Settings
        Settings of the simulation
        
        The dictionnary is formatted like this:
            
            dictionnary = {
                'variable_1': array of shape (nz,ny,nx),
                ...
                'variable_n': array of shape (nz,ny,nx)}
            
    """
    
    # write corresponding XDMF meta data file
    # for further info visit
    # https://pymotw.com/2/xml/etree/ElementTree/create.html
    # http://www.xdmf.org/index.php/XDMF_Model_and_Format
    rootxdmf = ET.Element(" ")
    rootxdmf.append(ET.Comment('DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []'))

    xdmf = ET.Element('Xdmf')
    xdmf.set("version", "2.0")
    domain=ET.SubElement(xdmf, "Domain")
    grid=ET.SubElement(domain, "Grid")

    topology=ET.SubElement(grid,'Topology')
    topology.set("TopologyType","3DRectMesh")
    topology.set("Dimensions",str(mesh.nz)+" "+str(settings.slice_lambda2)+" "+str(mesh.nx))

    grid.extend(topology)
    grid.set("Name", "mesh")
    grid.set("GridType", "Uniform")

    geometry=ET.SubElement(grid, "Geometry")
    geometry.set("GeometryType","VXVYVZ")

    # write mesh in xdmf file
    dataItemZc=ET.SubElement(geometry, "DataItem")
    dataItemZc.set("Dimensions",str(mesh.nx))
    dataItemZc.set("Name", "Xc")
    dataItemZc.set("NumberType", "Float")
    dataItemZc.set("Precision", "8")
    dataItemZc.set("Format", "HDF")
    dataItemZc.text = fname + '.h5'+":/Xc"

    dataItemYc=ET.SubElement(geometry, "DataItem")
    dataItemYc.set("Dimensions",str(settings.slice_lambda2))
    dataItemYc.set("Name", "Yc")
    dataItemYc.set("NumberType", "Float")
    dataItemYc.set("Precision", "8")
    dataItemYc.set("Format", "HDF")
    dataItemYc.text = fname + '.h5'+":/Yc"
    
    dataItemXc=ET.SubElement(geometry, "DataItem")
    dataItemXc.set("Dimensions",str(mesh.nz))
    dataItemXc.set("Name", "Zc")
    dataItemXc.set("NumberType", "Float")
    dataItemXc.set("Precision", "8")
    dataItemXc.set("Format", "HDF")
    dataItemXc.text = fname + '.h5'+":/Zc"

    # write 3D arrays in xdmf file
    for key in dictionnary.keys():

        attribute=ET.SubElement(grid, "Attribute")
        attribute.set("Name", str(key))
        attribute.set("AttributeType", "Scalar")
        attribute.set("Center","Node")
        dataItem=ET.SubElement(attribute, "DataItem")
        dataItem.set("Dimensions",str(mesh.nz)+" "+str(settings.slice_lambda2)+" "+str(mesh.nx))
        dataItem.set("NumberType", "Float")
        dataItem.set("Precision", "8")
        dataItem.set("Format", "HDF")
        dataItem.text = fname + '.h5'+":/" + str(key)

    # create corresponding file name by replacing file suffix
    with open(fname + '.xmf', "w+") as f:
        print(prettify(xdmf), file=f)
        print

    # add declaration, workaround to ET
    declaration='<!DOCTYPE Xdmf SYSTEM "xdmf.dtd" []>\n'
    f = open(fname + '.xmf', "r")
    contents = f.readlines()
    f.close()

    contents.insert(1, declaration)

    f = open(fname + '.xmf', "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()
    
    