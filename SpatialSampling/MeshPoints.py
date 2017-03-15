import numpy as np
import SpatialSampling
import Visualization
import scipy.ndimage as ndimage
from scipy import spatial

def mesh_dist(image, threshold, samplesize):
    newimage = np.zeros_like(image)
    # the follwing converts the zeros array newimage into a binary array
    newimage[image < threshold] = 0
    newimage[image >= threshold] = 1
    print "the binary image threshold is:  " + str(threshold) 
    # perform multiple image openings
    newimage = ndimage.binary_opening(newimage, structure=np.ones((3,3,3))).astype(np.int)
    newimage = ndimage.binary_dilation(newimage,iterations=samplesize).astype(newimage.dtype) # 40 is Sample_Density 
    # the following line consumes a lot of computation time@
    newimage2_dist = ndimage.distance_transform_edt(newimage)
    return newimage2_dist
        
def mesh_dist2(image, threshold): # no binary dilation in this function
    newimage = np.zeros_like(image)
    # the follwing converts the zeros array newimage into a binary array
    newimage[image < threshold] = 0
    newimage[image >= threshold] = 1
    print "the binary image threshold is:  " + str(threshold) 
    # perform multiple image openings
    newimage = ndimage.binary_opening(newimage, structure=np.ones((3,3,3))).astype(np.int)
    # the following line consumes a lot of computation time@
    newimage2_dist = ndimage.distance_transform_edt(newimage)
        
    return newimage2_dist
        
def sampling_def(image1, dist_img, threshold, Sample_Density):
    Samples2 = SpatialSampling.DiskSamples(image1, dist_img, threshold)
    Samples2.GenSamples(Sample_Density, 80000) # This distance will be used to calculate strains
    point_arr2 = Samples2.samples
    sample_point2 = np.transpose(point_arr2)
    Visualization.DataVisulization(image1, threshold).scatterplot(point_arr2)
    return sample_point2