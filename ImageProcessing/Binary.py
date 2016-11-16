"""

image transformation operations

"""

import numpy as np
from skimage.filters import threshold_otsu, threshold_adaptive

def toBinary(volume, threshold):
    """
    volume --- 3D volume
    threshold --- float value
    """
    if type(volume) != np.ndarray:
        
        raise TypeError('the input must be numpy array!')
        
    binary1 = np.zeros_like(volume) # make a copy of array
    
    binary1[ volume > threshold ] = 1
    
    return binary1
    
def sliceThreshold(volume, block_size = 5):
    """
    convert slice into binary using adaptive local ostu method
    
    volume --- 3D volume
    block_size --- int value
    
    """
    if type(volume) != np.ndarray:
        
        raise TypeError('the input must be numpy array!')

    x, y, z = volume.shape
    
    segImg = np.empty_like(volume)
    
    for i in range(z):
        
        binary_adaptive = threshold_adaptive(volume[:,:,i], block_size, offset=0)
        
        segImg[:,:,i] = binary_adaptive
    
    return segImg
