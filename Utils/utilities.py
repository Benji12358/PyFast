import numpy as np
import ffmpeg
import os
import h5py

def print_infos(arrayToShow, txt):

    print('Min : ', np.min(arrayToShow))
    print('Max : ', np.max(arrayToShow))
    print('########################################')

def find_nearest(array, value):
    """Given a mesh array, return the closest index corresponding to a given position. 
    
    Parameters
    ----------
    array : 3D numpy array with shape (n)
        Mesh array of any direction
    value : Real
        Position wanted
    
    Returns
    -------
    Integer
        The closest index corresponding to the value in the array
            
    """
    
    idx = (np.abs(array - value)).argmin()
    return idx


def read_h5_file(path, field_to_read, streamwise_direction):
    """Given a mesh array, return the closest index corresponding to a given position. 
    
    Parameters
    ----------
    array : 3D numpy array with shape (n)
        Mesh array of any direction
    value : Real
        Position wanted
    
    Returns
    -------
    Integer
        The closest index corresponding to the value in the array
            
    """
    
    f = h5py.File(path + '/' + field_to_read + '.h5', 'r')
    field = np.array(f[field_to_read])[:-1,:-1,:-1]
    f.close()

    if (streamwise_direction==3):
        field = np.transpose(field)

    return field


def getPathPicture(path, i):
    
    tmp=str(i)
    if i<10**(index-1):
        for j in range(len(str(i))-1,index-1):
            tmp = str(0) + tmp
    newPath = path + tmp + '.png'
    return newPath
    

def useffmpeg(picturesPath, videoPath):
    
    try:
        (
            ffmpeg
            .input(picturesPath, framerate=int((number_iteration-85)/10))
            .output(videoPath, acodec='libx264')
            .run(capture_stdout=True, capture_stderr=True)
        )
    except ffmpeg.Error as e:
        print('stdout:', e.stdout.decode('utf8'))
        print('stderr:', e.stderr.decode('utf8'))
        a=1

def createPostProcessingDirectory(path):
    """Create a directory of a given path, if it does not exist.  
    
    Parameters
    ----------
    path : String
        Path of the directory to create
            
    """
    
    if (not os.path.exists(path)):
        os.mkdir(path)