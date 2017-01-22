"""
This class/module saves the data into a folder or load a data from the folder.
"""
import numpy as np
import scipy.io as sio

def saveToDirectory(directoryFileName,nameInMat,par1):
    if type(directoryFileName) == str:
        sio.savemat(directoryFileName, {'par1':par1})
    else:
        print "Please check the format of input parameters!"
        
