"""
This is a tutorial for calibrating the corner solder balls
"""

import InputOutput  # JW
import SpatialSampling # JW
import Visualization # JW
#import CorrelationFunction #JW

import scipy.io as sio # pub
import scipy.ndimage as ndimage # pub
from mayavi import mlab # pub
import mayavi # pub
from scipy import spatial # pub
import numpy as np # pub
from scipy import spatial # pub
from scipy import linalg  # pub


#####  import the data from the local directory ###
PathDicom1 = "/Users/junchaowei/Desktop/kyber/Tutorial/Corner_Solder_Cal"
PathDicom2 = "/Users/junchaowei/Desktop/kyber/Tutorial/Corner_Solder_Cal"

dicom_volume1 = InputOutput.ReadDicom(PathDicom1)
dicom_volume2 = InputOutput.ReadDicom(PathDicom2)
            
dicom_volume1.loadfiles()
dicom_volume2.loadfiles()
    
volume1 = dicom_volume1.DicArray # raw

## filtering 
image1 = ndimage.filters.gaussian_filter(dicom_volume1.DicArray,3) # smoothed

threshold1 = 12000

def mesh_dist(image, threshold):
    """
    this function creates a distance image, basing on the Euclidean Distance from 
    the surface, identified by threshold value.
    input:
        image --- (N,d) np.array
        threshold --- a float value
    output:
        newimage2_dist --- (N,d) np.array
    """
    newimage = np.zeros_like(image)
    # the follwing converts the zeros array newimage into a binary array
    newimage[image < threshold] = 0
    newimage[image >= threshold] = 1
    print "the binary image threshold is:  " + str(threshold) 
    # perform multiple image openings
    ndimage.binary_opening(newimage, structure=np.ones((3,3,3))).astype(np.int)
    # the following line consumes a lot of computation time@
    newimage2_dist = ndimage.distance_transform_edt(newimage)
    return newimage2_dist
        
def sampling_def(image1, dist_img, threshold):
    """
    this is the sampling function
    
    input:
        image1 --- (N,d) np.array
        dist_img --- (N,d) np.array
        threshold --- a float value
    output:
        sample_point2 ---(N,3) np.array
    """
    Samples2 = SpatialSampling.DiskSamples(image1, dist_img, threshold)
    Samples2.GenSamples(3, 80000) # this number controls the density
    point_arr2 = Samples2.samples
    sample_point2 = np.transpose(point_arr2)
    Visualization.DataVisulization(image1, threshold).scatterplot(point_arr2)
    return sample_point2

dist_img1 = mesh_dist(image1, threshold1)
point_arr1 = sampling_def(image1, dist_img1, threshold1)

vl = Visualization.DataVisulization(image1, threshold1)
vl.contour3d()



 