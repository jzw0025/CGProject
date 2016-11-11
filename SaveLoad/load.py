import numpy as np
import scipy.io as sio

def loadToDirectory(directoryFileName):
    mat_contents1 = sio.loadmat(directoryFileName)
    data = mat_contents1['par1'] # default saved data name.
    return data