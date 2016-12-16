"""
this cleans the image, and separates the big image into pieces

"""

import scipy.io as sio
import matplotlib.pyplot as plt
import ImageProcessing
import numpy as np
import Visualization
from scipy import ndimage

volume000 = sio.loadmat('/Users/junchaowei/Desktop/SpaceRegistration_000_125/volume000.mat') # read the file     
im1 = volume000['par1']#[280:500,90:350,:] #[350:450,150:250,:]                                                                                                                                                                                                                                                                                                       
volume125 = sio.loadmat('/Users/junchaowei/Desktop/SpaceRegistration_000_125/volume125_regi.mat')# big region three point alignment 
im2 = volume125['par1']#[280:500,90:350,:]

fig = plt.figure()
plt.imshow(im1[:,:,70],cmap='Greys_r')  
ax = fig.add_subplot(111)
ax.set_title('click to build line segments')
line, = ax.plot([0], [0])  # empty line
linebuilder = ImageProcessing.LineBuilder(line)  
plt.show()

x_cor = linebuilder.xs
y_cor = linebuilder.ys

points = np.array([x_cor, y_cor]).T
points = np.delete(points, (0), axis=0)

imx, imy = im1[:,:,60].shape
mx, my = np.meshgrid(np.linspace(0,imx-1,imx), np.linspace(0,imy-1,imy))
# meshgrid change the orientation of x,y: x becomes the horizontal axis, y becomes the vertical
# this change would affect the return statement of reshape.
Cord_xy = np.vstack([mx.flatten(), my.flatten()]).T # column major vector

#### run previous code to choose the points
Inhull_point = ImageProcessing.in_hull(Cord_xy, points)
Index_logic = Inhull_point.reshape(imx, imy)
plt.imshow(Index_logic)

### take the region 1 out from big image
region1_ref = np.empty_like(im1)
for i in range(im1.shape[2]):
    region1_ref[:,:,i] = im1[:,:,i]*Index_logic
    
### take the region 2 out from big image
region1_def = np.empty_like(im2)
for i in range(im2.shape[2]):
    region1_def[:,:,i] = im2[:,:,i]*Index_logic
    
region1_def_hist = ImageProcessing.hist_match(region1_def, region1_ref)
    
plot = Visualization.DataVisulization(ndimage.gaussian_filter(region1_ref,5), 90)
plot.contour3d()
plot = Visualization.DataVisulization(ndimage.gaussian_filter(region1_def_hist,5), 90)
plot.contour3d()

sio.savemat("/Users/junchaowei/Desktop/Region1/region1_ref.mat", {"par1":region1_ref})
sio.savemat("/Users/junchaowei/Desktop/Region1/region1_def.mat", {"par1":region1_def_hist})