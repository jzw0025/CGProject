"""
this is a developing technique for image registration

"""
from mayavi import mlab
from scipy import ndimage
import numpy as np
import array
import scipy.io as sio
import Visualization
import scipy.ndimage as ndimage
import ImageProcessing
from scipy import linalg
import SaveLoad
import matplotlib.pyplot as plt
from scipy import interpolate

mlab.figure(bgcolor=(0.5,0.5,0.5), size=(400, 400)) # set the background color (r,g,b) form [0.0,1.0]

volume000 = sio.loadmat('F:\For JC\practice run/volume000.mat') # read the file     
im1 = volume000['par1']#[280:500,90:350,:] #[350:450,150:250,:]                                                                                                                                                                

volume125 = sio.loadmat('F:\For JC\practice run/volume125.mat')# big region three point alignment 
im2 = volume125['par1']#[280:500,90:350,:]

im2_matched = ImageProcessing.hist_match(im2, im1) 

plot = Visualization.DataVisulization(ndimage.gaussian_filter(im1,5), 70)
plot.contour3d()
plot2 = Visualization.DataVisulization(ndimage.gaussian_filter(im2_matched,5), 70)
plot2.contour3d()


plt.figure(1)
plt.imshow(im1[:,:,60], cmap = "Greys_r")
plt.show()

sy,sx = im1[:,:,0].shape
registration_check = np.zeros((sy,sx,3))
registration_check[:,:,0] = 255-im1[:,:,60]
registration_check[:,:,1] = 255-im2[:,:,60]#.astype(int) # the float value must be converted into int color channel
registration_check[:,:,2] = np.zeros((sy,sx))
plt.figure(99)
plt.imshow(registration_check)
plt.show()

nodes, elements = plot.surf_points, plot.surf_elements
nodes2, elements2 = plot2.surf_points, plot2.surf_elements

whole_nodes = np.concatenate((nodes,nodes2), axis=0)
sx, sy = nodes.shape
elements3 = elements2 + sx
whole_elements = np.concatenate((elements, elements3), axis=0)

#### help function #####
def reshapePoint(point1, point2): # reshape the points array
    xpos = np.ndarray((2)) 
    ypos = np.ndarray((2)) 
    zpos = np.ndarray((2)) 
    xpos[0] = point1[0]
    xpos[1] = point2[0]
    ypos[0] = point1[1]
    ypos[1] = point2[1]
    zpos[0] = point1[2]
    zpos[1] = point2[2]
    color = 0.5*np.ones(xpos.shape) 
    return xpos, ypos, zpos, color

# A first plot in 3D
fig = mlab.figure(2)
mlab.clf()
mesh = mlab.triangular_mesh(whole_nodes[:,0], whole_nodes[:,1], whole_nodes[:,2], whole_elements)
cursor3d = mlab.points3d(0., 0., 0., mode='axes',
                                color=(0, 0, 0),
                                scale_factor=0.5)
                                
mlab.title('Click on the LC Region')

## A second plot, flat projection image, this is under development
#fig2d = mlab.figure(2)
#mlab.clf()
#im = mlab.imshow(s)
#cursor = mlab.points3d(0, 0, 0, mode='2dthick_cross',
#                                color=(0, 0, 0),
#                                scale_factor=10)
#mlab.view(90, 0)

start_points = []
end_points = []
counter = 1
################################################################################
def picker_callback(picker_obj):
    global counter, start_points, end_points
    counter += 1
    picked = picker_obj.actors
    print "counter call:" + str(counter)
    if mesh.actor.actor._vtk_obj in [o._vtk_obj for o in picked]:
        # m.mlab_source.points is the points array underlying the vtk
        # dataset. GetPointId return the index in this array.
        print picker_obj.point_id
        x_ = np.lib.index_tricks.unravel_index(picker_obj.point_id,whole_nodes[:,0].shape)                                                             
        print whole_nodes[:,0][x_], whole_nodes[:,1][x_], whole_nodes[:,2][x_]
        cursor3d.mlab_source.reset(x=whole_nodes[:,0][x_],
                                   y=whole_nodes[:,1][x_],
                                   z=whole_nodes[:,2][x_])                          
fig.on_mouse_pick(picker_callback)
mlab.show()

def fixed_SVD(target, source):
    """
    this module calculates the fixed SVD for the paired input points
    
    """
    translation = np.zeros([1,3])
    vector = target.mean(0) - source.mean(0) # first dimension mean
    source = source + vector
    ##########  rotation  ############
    matrix1 = np.matrix(source)
    matrix2 = np.matrix(target)
    # building the corresponding matrix
    resultMatrix = matrix1.T.dot(matrix2)
    # perform the SVD algorithm
    U, s, Vh = linalg.svd(resultMatrix)
    Rotation = Vh.T.dot(U.T)
    source2 = Rotation.dot(np.matrix(source).T).T
    #print newsource == newsource2
    source = np.array(source2) # updating the source points
    translation = translation + vector # translation is zero initially, this does not need to add iteratively, but keep the original formula
            
    return source, Rotation, translation 
    



target = [587.64642334 ,327.0 ,12.0
        80.0, 169.0, 26.866645813
        99.0, 403.0, 26.5643672943]]

source = [
          [172.88319397, 55.0, 65.0],
          [86.0, 106.616584778, 82.0],
          [22.2802658081, 116.0, 38.0]]
                 
new = np.concatenate((np.array(target).T,np.array(source).T),axis=1)
Visualization.DataVisulization(im1, 1).scatterplot(new)

newPoints, rotation2, translation2 = fixed_SVD(np.array(target), np.array(source))

new2 = np.concatenate((np.array(target).T,np.array(newPoints).T),axis=1)
Visualization.DataVisulization(im1, 1).scatterplot(new2)

translation3 =  ((np.array(source).mean(0)).dot(rotation2)-np.array(target).mean(0)).dot(linalg.inv(rotation2))
newim3 = ndimage.interpolation.affine_transform(im2_matched,rotation2,order=3,offset=translation3,cval=0.0) # rigid affine registration 

##### Affine Model Displacement #########
affineMatrix = np.zeros((4,4))
affineMatrix[0:3,0:3] = rotation2
affineMatrix[0:3,3] = translation3
affineMatrix[3,3] = 1.0
ix, iy, iz = im1.shape
mapX = np.empty_like(im1)
mapY = np.empty_like(im1)
mapZ = np.empty_like(im1)
displacementX = np.empty_like(im1)
displacementY = np.empty_like(im1)
displacementZ = np.empty_like(im1)
for i in range(ix):
    for j in range(iy):
        for k in range(iz):
            # coordinator is (i,j,k,1)
            cord = np.array([i,j,k,1])
            new_cord = affineMatrix.dot(cord)
            mapX[i,j,k] = new_cord[0]
            mapY[i,j,k] = new_cord[1]
            mapZ[i,j,k] = new_cord[2]
            displacementX[i,j,k] = new_cord[0]-i
            displacementY[i,j,k] = new_cord[1]-j
            displacementZ[i,j,k] = new_cord[2]-k
            
newim3_map = ndimage.map_coordinates(im2_matched, [mapX, mapY, mapZ], order=1)

nodeDisplacementX = np.empty(nodes.shape[0])        
nodeDisplacementY = np.empty(nodes.shape[0]) 
nodeDisplacementZ = np.empty(nodes.shape[0]) 
for inode in range(nodes.shape[0]):
    cord = np.array([nodes[inode][0],nodes[inode][1],nodes[inode][2],1])
    new_cord = affineMatrix.dot(cord)
    nodeDisplacementX[inode] = new_cord[0]-nodes[inode][0]
    nodeDisplacementY[inode] = new_cord[1]-nodes[inode][1]
    nodeDisplacementZ[inode] = new_cord[2]-nodes[inode][2]

surf = plot.surf_points
inter_z = interpolate.LinearNDInterpolator(nodes, nodeDisplacementY, fill_value=0.0, rescale=True)        
extrapo_z = inter_z.__call__(surf)
plot.plot3dpoint(extrapo_z)
mlab.colorbar(object=None, title="test", orientation="vertical",nb_labels=20, nv_colors=None, label_fmt=None)

#####
plot3 = Visualization.DataVisulization(ndimage.gaussian_filter(newim3_map,5), 110)
plot3.contour3d()
plot4 = Visualization.DataVisulization(ndimage.gaussian_filter(newim3,5), 110)
plot4.contour3d()


save = SaveLoad.saveDataBase("C:\Users\HT\Desktop/test run/") 
datadict={}
datadict['volume000'] = nodeDisplacementX
datadict['volume125'] = nodeDisplacementX
datadict['volume125_regi'] = nodeDisplacementX
datadict['target_point'] = im1
datadict['source_point'] = im2_matched
datadict['source_point'] = newim3
datadict['target_point'] = target
datadict['source_point'] = source
save.save(datadict)


###### multi-channel color images ####
sy,sx = im1[:,:,0].shape
registration_check = np.zeros((sy,sx,3))
registration_check[:,:,0] = 255-im1[:,:,70]
registration_check[:,:,1] = 255-newim3[:,:,70].astype(int) # the float value must be converted into int color channel
registration_check[:,:,2] = np.zeros((sy,sx))
plt.figure(99)
plt.imshow(registration_check)
plt.show()

################ looping the figures and saving  #############
ix,iy,iz = im1.shape 
for i in range(iz):
    registration_check = np.zeros((sy,sx,3))
    registration_check[:,:,0] = 255-im1[:,:,i]
    registration_check[:,:,1] = 255-newim3[:,:,i].astype(int) # the float value must be converted into int color channel
    registration_check[:,:,2] = np.zeros((sy,sx))
    plt.figure(i)
    plt.imshow(registration_check)
    plt.show()
    plt.savefig('C:\Users\HT\Desktop/test run/stack/' + str(i) + '.png')
    plt.close()

