"""
this is a file for combining results from different regions of the same image

"""
import scipy.io as sio
import numpy as np


volume = sio.loadmat("C:/Users/HT/Desktop/test run/combine regions/Region1/mask.mat")
mask1 = volume["par1"]
volume = sio.loadmat("C:/Users/HT/Desktop/test run/combine regions/Region1/displacementZ.mat")
region1_displacementZ = volume["par1"]
volume = sio.loadmat("C:/Users/HT/Desktop/test run/combine regions/Region1/region1_ref.mat")
image1 = volume["par1"]


volume = sio.loadmat("C:/Users/HT/Desktop/test run/combine regions/Region2/mask.mat")
mask2 = volume["par1"]
volume = sio.loadmat("C:/Users/HT/Desktop/test run/combine regions/Region2/displacementZ.mat")
region2_displacementZ = volume["par1"]

print mask1.shape
print mask2.shape
print region1_displacementZ.shape
print region2_displacementZ.shape

sx,sy,sz = region1_displacementZ.shape

for i in range(sz):
    print "process slice:" + str(i)
    region1_displacementZ[:,:,i] = region1_displacementZ[:,:,i]*mask1
    region2_displacementZ[:,:,i] = region2_displacementZ[:,:,i]*mask2
    
    
overall_displacementZ = np.empty_like(region1_displacementZ)

overall_displacementZ = overall_displacementZ + region1_displacementZ
overall_displacementZ = overall_displacementZ + region2_displacementZ


import matplotlib.pyplot as plt

plt.figure(1)
plt.imshow(overall_displacementZ[:,:,90])
plt.show()

from scipy import ndimage
import Visualization
from scipy import interpolate

plot = Visualization.DataVisulization(ndimage.gaussian_filter(image1,5), 110)
plot.contour3d()

surf = plot.surf_points
inter_z = interpolate.LinearNDInterpolator(nodes, nodeDisplacementY, fill_value=0.0, rescale=True)        
extrapo_z = inter_z.__call__(surf)
plot.plot3dpoint(extrapo_z)
mlab.colorbar(object=None, title="test", orientation="vertical",nb_labels=20, nv_colors=None, label_fmt=None)

