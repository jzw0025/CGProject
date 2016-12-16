"""
proces regions
"""
   
import scipy.io as sio
import ImageProcessing
import Visualization
from scipy import ndimage
import matplotlib.pyplot as plt

volume000 = sio.loadmat('C:/Users/HT/Desktop/test run/combine regions/Region1/region1_ref.mat') # read the file     
im1 = volume000['par1']                                                                                                                                                                                                                                                                                           
volume125 = sio.loadmat('C:/Users/HT/Desktop/test run/combine regions/Region1/region1_def.mat')# big region three point alignment 
im2 = volume125['par1']
volume = sio.loadmat('C:/Users/HT/Desktop/test run/combine regions/Region1/mask.mat')# big region three point alignment 
mask = volume['par1']

plt.imshow(mask)
plt.show()

#plot = Visualization.DataVisulization(ndimage.gaussian_filter(im1,5), 70)
#plot.contour3d()
#plot2 = Visualization.DataVisulization(ndimage.gaussian_filter(im2,5), 70)
#plot2.contour3d()

LC = ImageProcessing.regionLC(im1,im2)

LC.selectPoints(5,110)

LC.addPairPoints([450.0,151.621856689,68.0],[451.0,143.87109375,66.0])

LC.addPairPoints([340.0,126.921035767,70.0],[344.480072021,121.0,65.0])

LC.addPairPoints([362.0,206.564971924,79.0],[367.857025146,199.0,75.0])


LC.computeAffine()

sio.savemat("C:/Users/HT/Desktop/test run/combine regions/Region1/target_list", {"par1":LC.target_list})
sio.savemat("C:/Users/HT/Desktop/test run/combine regions/Region1/source_list", {"par1":LC.source_list})

sio.savemat("C:/Users/HT/Desktop/test run/combine regions/Region1/region2_regi.mat", {"par1":LC.image3})

sio.savemat("C:/Users/HT/Desktop/test run/combine regions/Region1/displacementX", {"par1":LC.displacementX})
sio.savemat("C:/Users/HT/Desktop/test run/combine regions/Region1/displacementY", {"par1":LC.displacementY})
sio.savemat("C:/Users/HT/Desktop/test run/combine regions/Region1/displacementZ", {"par1":LC.displacementZ})

