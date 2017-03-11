"""
region class
"""
import Visualization
from scipy import ndimage
import numpy as np
from mayavi import mlab
from scipy import linalg

class regionLC():
    
    def __init__(self, regionImage1, regionImage2):
        self.image1 = regionImage1
        self.image2 = regionImage2
        self.image3 = []
        self.target_list = []
        self.source_list = []
        self.displacementX = []
        self.displacementY = []
        self.displacementZ = []
        
    def selectPoints(self, visu_factor, iso):
        
        plot = Visualization.DataVisulization(ndimage.gaussian_filter(self.image1,visu_factor), iso)
        plot.contour3d()
        
        plot2 = Visualization.DataVisulization(ndimage.gaussian_filter(self.image2,visu_factor), iso)
        plot2.contour3d()
        
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
        
        ################################################################################
        def picker_callback(picker_obj):
            print "selecting a new point"
            picked = picker_obj.actors
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
        
    def addPairPoints(self,target, source):
        if len(target)<3:
            print "warning: at least three points are required to compute the transformation!"
        
        print "current length of target points list: " + str(len(self.target_list))
        
        self.target_list.append(target)
        self.source_list.append(source)
            
    def computeAffine(self):
        
        ###  this is help function ###
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
            
        new = np.concatenate((np.array(self.target_list).T,np.array(self.source_list).T),axis=1)
        Visualization.DataVisulization(self.image1, 1).scatterplot(new)
        
        newPoints, rotation2, translation2 = fixed_SVD(np.array(self.target_list), np.array(self.source_list))
        
        new2 = np.concatenate((np.array(self.target_list).T,np.array(newPoints).T),axis=1) # visualize transformed source point and target point
        
        Visualization.DataVisulization(self.image1, 1).scatterplot(new2)
        
        translation3 =  ((np.array(self.source_list).mean(0)).dot(rotation2)-np.array(self.target_list).mean(0)).dot(linalg.inv(rotation2))
        self.image3 = ndimage.interpolation.affine_transform(self.image2,rotation2,order=3,offset=translation3,cval=0.0) # rigid affine registration 
        
        ##### Affine Model Displacement #########
        affineMatrix = np.zeros((4,4))
        affineMatrix[0:3,0:3] = rotation2
        affineMatrix[0:3,3] = translation3
        affineMatrix[3,3] = 1.0
        
        ix, iy, iz = self.image1.shape
        
        mapX = np.empty_like(self.image1)
        mapY = np.empty_like(self.image1)
        mapZ = np.empty_like(self.image1)
        
        self.displacementX = np.empty_like(self.image1)
        self.displacementY = np.empty_like(self.image1)
        self.displacementZ = np.empty_like(self.image1)
        
        for i in range(ix):
            for j in range(iy):
                for k in range(iz):
                    # coordinator is (i,j,k,1)
                    cord = np.array([i,j,k,1])
                    new_cord = affineMatrix.dot(cord)
                    mapX[i,j,k] = new_cord[0]
                    mapY[i,j,k] = new_cord[1]
                    mapZ[i,j,k] = new_cord[2]
                    self.displacementX[i,j,k] = new_cord[0]-i
                    self.displacementY[i,j,k] = new_cord[1]-j
                    self.displacementZ[i,j,k] = new_cord[2]-k
                    
    def getDisplacementX(self):
        return self.displacementX
        
    def getDisplacementY(self):
        return self.displacementY
        
    def getDisplacementZ(self):
        return self.displacementZ
    
    def getTargetList(self):
        return self.target_list
    
    def getSourceList(self):
        return self.source_list
                    
if __name__ == "__main__":
    print "this is a unit test file"
    import scipy.io as sio
    volume000 = sio.loadmat('/Users/junchaowei/Desktop/Pack12122016/clean_000.mat') # read the file     
    im1 = volume000['par1'][200:300,200:300,50:100]                                                                                                                                                                                                                                                                                                    
    volume125 = sio.loadmat('/Users/junchaowei/Desktop/Pack12122016/clean487_reg1.mat')# big region three point alignment 
    im2 = volume125['par1'][200:300,200:300,50:100]
    
    LC =  regionLC(im1,im2)  
    
    LC.selectPoints(5,im1.mean())
    
    #LC.addPairPoints([95.0,120.403953552,100.0],[102.0,113.3279953,90.0])
    #LC.addPairPoints([169.642150879,62.0,70.0],[171.0,54.0,66.319732666])
    #LC.addPairPoints([42.3362579346,127.0,89.0],[46.5657348633,120.0,83.0])
    
    for i in range(len(pairs)):
        LC.addPairPoints(pairs[i][0],pairs[i][1])
        
    LC.computeAffine()
          
    sio.savemat("/Users/junchaowei/Desktop/Region2/target_list", {"par1":LC.target_list})
    sio.savemat("/Users/junchaowei/Desktop/Region2/source_list", {"par1":LC.source_list})
    
    sio.savemat("/Users/junchaowei/Desktop/Region2/image1", {"par1":LC.image1})
    sio.savemat("/Users/junchaowei/Desktop/Region2/image2", {"par1":LC.image2})
    sio.savemat("/Users/junchaowei/Desktop/Region2/image2_regi", {"par1":LC.image3})

    sio.savemat("/Users/junchaowei/Desktop/Region2/displacementX", {"par1":LC.displacementX})
    sio.savemat("/Users/junchaowei/Desktop/Region2/displacementY", {"par1":LC.displacementY})
    sio.savemat("/Users/junchaowei/Desktop/Region2/displacementZ", {"par1":LC.displacementZ})

               