"""
convert the volume into slice images into a fold

"""

import numpy as np
import matplotlib.pyplot as plt

def ConvertImageSlice(volume, address=None):
    """
    volume --- ndarray 
    address --- folder address
    
    """
    if not address:
        print "no address is imported: please make a fold on desktop named: image"
        address = '/Users/junchaowei/Desktop/image/'
    
    for i in range(volume.shape[2]):
        print "converting slice:" + str(i)
        plt.imshow(volume[:,:,i], interpolation='nearest', cmap='Greys_r')
        plt.savefig(address+ str(i) +'.tif')
        
    print "finished converting image slices"
    

if __name__ == "__main__":
    print "this is a unit test!"
    
    address = "/Users/junchaowei/Desktop/deformed000_cleanLC/"
    ConvertImageSlice(dynamicRangeImage, address)
    
    address2 = "/Users/junchaowei/Desktop/deformed487_cleanLC/"
    ConvertImageSlice(dynamicRangeImage2, address2)