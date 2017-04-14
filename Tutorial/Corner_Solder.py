"""
This is a tutorial for calibrating the corner solder balls
"""

import InputOutput  # JW
import SpatialSampling # JW
import Visualization # JW
import CorrelationFunction #JW

import scipy.io as sio # pub
import scipy.ndimage as ndimage # pub
from mayavi import mlab # pub
import mayavi # pub
from scipy import spatial # pub
import numpy as np # pub
from scipy import spatial # pub
from scipy import linalg  # pub
from scipy import interpolate


#####  import the data from the local directory ###
PathDicom1 = "/Corner_Solder_Cal"
PathDicom2 = "/Corner_Solder_Cal"

dicom_volume1 = InputOutput.ReadDicom(PathDicom1)
dicom_volume2 = InputOutput.ReadDicom(PathDicom2)
            
dicom_volume1.loadfiles()
dicom_volume2.loadfiles()
    
volume1 = dicom_volume1.DicArray # raw
volume2 = dicom_volume2.DicArray # raw

## filtering 
image1 = ndimage.filters.gaussian_filter(dicom_volume1.DicArray,3) # smoothed
image2 = ndimage.filters.gaussian_filter(dicom_volume2.DicArray,3) # smoothed

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

dist_img2 = mesh_dist(image2, threshold1)
point_arr2 = sampling_def(image2, dist_img2, threshold1)

vl = Visualization.DataVisulization(image1, threshold1)
vl.contour3d()

Sample_Density = 3
mesh2 = SpatialSampling.MeshDelaunay(point_arr1.T) # use the original point_arr1 
mesh2.mesh()
mesh2.alpha_shape(Sample_Density)
mesh2.mesh_triangle() # find unique triangles
    
for dummy_i in range(1):   
    Mesh_optimizer = SpatialSampling.MeshOptimizer(point_arr1.T,mesh2.ntri)
    Mesh_optimizer.edge_connector_smooth()
  
NodeList = []
subSize = 21
subSearch = 15
dictions = CorrelationFunction.calculateSubset(image1, image2, mesh2.points, subSize, subSearch)
sub1 = dictions[0] 
nodesDiction={}  

for i in range(mesh2.points.shape[0]):
        print i
        tempNode = CorrelationFunction.Node([mesh2.points[i][0], mesh2.points[i][1], mesh2.points[i][2]], subSize, subSearch, i)
        tempNode.addSubset(dictions[0].get(i), dictions[1].get(i), dictions[4].get(i)) 
        # the output index dictions[index]: 0-subset1, 1-subset2, 2-correlation value, 3-displacements, 4-psubset2
        # add the correlation value to each node
        tempNode.addCorrelation(dictions[2].get(i))
        # add initial search vector to each node
        tempNode.addVector(dictions[3].get(i)[0], dictions[3].get(i)[1], dictions[3].get(i)[2]) # from dictionary adding three displacements (z,x,y) value into a node instance
        ## calculating connecting edge point to a point
        connected_points = CorrelationFunction.ConnectMesh(mesh2.ntri,i)
        # add the connectivity to the node
        ##tempNode.addConnectivity(connected_points)
        # append the node to the list
        NodeList.append(tempNode)
        # put node in the dictionary for later search
        nodesDiction[i] = tempNode
        # sort the node according to the correlation value
     
index = 0
outputColor = []
size = str(len(NodeList))
    
def ICP_KLT(point1, point2, iterations):
        def SVD(target, source):
            """
            the input data:
            target --- N-D array e.g. shape: (1000, 3)
            source --- N-D array e.g. shape: (1000, 3) 
            the return data:
            
            """
            targetTree = spatial.KDTree(target) # create K-D tree for the input target 
            ######### the following code goes iterately ##########
            translation = np.zeros([1,3])
            ##########   translation  ########
            vector = target.mean(0) - source.mean(0)
            source = source + vector
            ##########  rotation  ############
            sourcePair_list = []
            targetPair_list = []       
            for i in range(source.shape[0]):
                # find the closest point in the tree-data
                results = targetTree.query(source[i])
                sourcePair_list.append(list(source[i]))
                targetPair_list.append(list(target[results[1]]))           
            matrix1 = np.matrix(sourcePair_list)
            matrix2 = np.matrix(targetPair_list)
            # building the corresponding matrix
            resultMatrix = matrix1.T.dot(matrix2)
            # perform the SVD algorithm
            U, s, Vh = linalg.svd(resultMatrix)
            Rotation = Vh.T.dot(U.T)
            source2 = Rotation.dot(np.matrix(source).T).T
            #print newsource == newsource2
            source = np.array(source2) # updating the source points
            translation = translation + vector
            return source, Rotation, translation 
        
        change_points = point2
        translation2 = [0.0,0.0,0.0]
        for i in range(iterations): # this SVD registration is followed by image registration at each step.
            old_change_points = change_points
            change_points, R, V = SVD(point1, change_points)
            print "the sum sqare for the point cloud is: " + str(np.sqrt(sum(sum((change_points-old_change_points)**2))))
            rotation2 = R
            #translation = self.mean1_vector.T[0] - np.dot(np.dot(self.eig_vec_cov1,self.eig_vec_cov2.T),self.mean2_vector.T[0])
            #translation2 = ((point2.mean(0)).dot(rotation2)-point1.mean(0)).dot(linalg.inv(rotation2))   
            translation2 = translation2 + V          
        return translation2 
        
def subPoint(point_arr, center, radius):  
               
        x1 = center[0] - radius
        x2 = center[0] + radius 
        
        y1 = center[1] - radius
        y2 = center[1] + radius   
        
        z1 = center[2] - radius
        z2 = center[2] + radius
        
        index_x = np.logical_and(point_arr[:,0]>=x1, point_arr[:,0]<=x2)
        index_y = np.logical_and(point_arr[:,1]>=y1, point_arr[:,1]<=y2)
        index_z = np.logical_and(point_arr[:,2]>=z1, point_arr[:,2]<=z2)    
        index = np.logical_and(index_x, index_y, index_z)
        
        return point_arr[index]

for eachNote in NodeList:
        print "the computation process is: "+str(NodeList.index(eachNote))+ " / "+ size
        print "the correlation value of node is: "+ str(eachNote.correlationValue)
        
        #if i == 100:
        #    break
        #i+=1
        
        if eachNote.initialGuess == False: # if the node has not been calculated.
            #subsets = eachNote.getKLTsubset()
            
            old_loc = eachNote.getOldCoordinator()
            new_loc = eachNote.getNewCoordinator()     
             
            image1_coordinator = np.ceil(np.array([old_loc[0], old_loc[1], old_loc[2]]))
            image2_coordinator = np.ceil(np.array([new_loc[0], new_loc[1], new_loc[2]]))  
                      
            ref_points = subPoint(point_arr1, image1_coordinator, 5*(subSize-1)/2)
            def_points = subPoint(point_arr2, image2_coordinator, 5*(subSize-1)/2) 
            
            #ref_points[:,0] = ref_points[:,0] + eachNote._initZ
            #ref_points[:,1] = ref_points[:,1] + eachNote._initX
            #ref_points[:,2] = ref_points[:,2] + eachNote._initY
                  
            translation = ICP_KLT(def_points, ref_points, 4)
            eachNote.opt_x = translation[0][0]
            eachNote.opt_y = translation[0][1]
            eachNote.opt_z = translation[0][2]
            eachNote.initialGuess = True      
            #print eachNote.getConnectivity()   
            #break

dispx = np.zeros(len(nodesDiction))
dispxx = np.zeros(len(nodesDiction))
dispy = np.zeros(len(nodesDiction))
dispyy = np.zeros(len(nodesDiction))
dispz = np.zeros(len(nodesDiction))
dispzz = np.zeros(len(nodesDiction))
nodes_dis = []
    
for icolorNode in nodesDiction.keys():
        dispx[icolorNode] = nodesDiction.get(icolorNode).opt_x
        dispxx[icolorNode] = nodesDiction.get(icolorNode)._initZ
        dispy[icolorNode] = nodesDiction.get(icolorNode).opt_y
        dispyy[icolorNode] = nodesDiction.get(icolorNode)._initX
        dispz[icolorNode] = nodesDiction.get(icolorNode).opt_z
        dispzz[icolorNode] = nodesDiction.get(icolorNode)._initY
        nodes_dis.append(nodesDiction.get(icolorNode).getOldCoordinator())
     
    #dispz = np.zeros(mesh2.points.T.shape[1])
colors = (dispz)*0.0495
    #colors = np.zeros(mesh2.points.T.shape[1])
    #sio.savemat('/Users/junchaowei/Desktop/0_90_z.mat', {'p':colors})
    #sio.savemat('/Users/junchaowei/Desktop/0_90_element.mat',{'p':mesh2.ntri})
    
@mlab.show
def main():
    mesh2.view(mesh2.polydata(mesh2.points.T,colors))  
main()
mayavi.mlab.colorbar(object=None, title="U(Z)(mm)", orientation="vertical", nb_labels=20, nb_colors=None, label_fmt=None)
"colorbar reference: http://docs.enthought.com/mayavi/mayavi/auto/mlab_decorations.html"
    
plot = Visualization.DataVisulization(image1, threshold1+5000)
plot.contour3d()
surf = plot.surf_points
inter_z = interpolate.LinearNDInterpolator(mesh2.points, colors, fill_value=0.0, rescale=True)        
extrapo_z = inter_z.__call__(surf)
plot.plot3dpoint(extrapo_z)
