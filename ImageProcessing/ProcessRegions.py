"""
proces regions
"""
   
import scipy.io as sio
import ImageProcessing
import Visualization
from scipy import ndimage
import matplotlib.pyplot as plt

volume000 = sio.loadmat('/Users/junchaowei/Desktop/Region2/region1_ref.mat') # read the file     
im1 = volume000['par1']                                                                                                                                                                                                                                                                                           
volume125 = sio.loadmat('/Users/junchaowei/Desktop/Region2/region1_def.mat')# big region three point alignment 
im2 = volume125['par1']
volume = sio.loadmat('/Users/junchaowei/Desktop/Region2/mask.mat')# big region three point alignment 
mask = volume['par1']

plt.imshow(mask)
plt.show()

#plot = Visualization.DataVisulization(ndimage.gaussian_filter(im1,5), 70)
#plot.contour3d()
#plot2 = Visualization.DataVisulization(ndimage.gaussian_filter(im2,5), 70)
#plot2.contour3d()

LC = ImageProcessing.regionLC(im1,im2)

LC.selectPoints(5,70)

LC.addPairPoints([272.0,134.635803223,83.0],[262.0,130.0,79.8894424438])
LC.addPairPoints([271.0,166.0,87.6977920532],[256.0,163.0,78.2035446167])
LC.addPairPoints([261.0,194.0,92.272102356],[251.229248047, 189.0, 89.0])
LC.addPairPoints([294.0,113.0,73.3494644165],[283.0,104.68258667,80.0])

LC.computeAffine()

sio.savemat("/Users/junchaowei/Desktop/Region2/target_list", {"par1":LC.target_list})
sio.savemat("/Users/junchaowei/Desktop/Region2/source_list", {"par1":LC.source_list})

sio.savemat("/Users/junchaowei/Desktop/Region2/region2_regi.mat", {"par1":LC.image3})

sio.savemat("/Users/junchaowei/Desktop/Region2/displacementX", {"par1":LC.displacementX})
sio.savemat("/Users/junchaowei/Desktop/Region2/displacementY", {"par1":LC.displacementY})
sio.savemat("/Users/junchaowei/Desktop/Region2/displacementZ", {"par1":LC.displacementZ})


