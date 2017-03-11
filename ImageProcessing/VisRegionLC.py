import Visualization
import scipy.io as sio
import numpy as np

load = sio.loadmat('/Users/junchaowei/Desktop/Region2/image1.mat')
image1 = load['par1']
load = sio.loadmat('/Users/junchaowei/Desktop/Region2/image2.mat')
image2 = load['par1']
load = sio.loadmat('/Users/junchaowei/Desktop/Region2/image2_regi.mat')
image2_r = load['par1']

v1 = Visualization.DataVisulization(image1,image1.mean()) # initial reference image
v1.contour3d()
v1 = Visualization.DataVisulization(image2,image2.mean()) # initial reference image
v1.contour3d()
v1 = Visualization.DataVisulization(image2_r,image2_r.mean()) # initial reference image
v1.contour3d()

import matplotlib.pyplot as plt
sy,sx = image1[:,:,0].shape
registration_check = np.zeros((sy,sx,3))
registration_check[:,:,0] = image1[:,:,25]
registration_check[:,:,1] = image2[:,:,25]#.astype(int) # the float value must be converted into int color channel
registration_check[:,:,2] = np.zeros((sy,sx))
plt.figure(100)
plt.imshow(registration_check)
plt.show()