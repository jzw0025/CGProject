"""
this class is the create an critical node.
each node should have the best displacement, and it stores 12 initial search constants

it also contains initial displacement vectors.

"""
import numpy as np
import CorrelationFunction 

class Node:
    
    def __init__(self, location, subSize, subSearch, number):
        
        self.number = number # this stores the node numbering
        
        self._z = location[0] # those are three coordinators of node.
        self._x = location[1]
        self._y = location[2]
        
        self._initZ = None
        self._initX = None
        self._initY = None
                
        ######## subset coordinators ##########
        self.subset = subSize/2
        self.subSearch = subSearch
        
        ########  KLT initial values  #########
        self.initialGuess = False
        
        self._u = 0.0 #three displacement vectors
        self._v = 0.0
        self._w = 0.0
        
        self._ux = 0.0 # partial derivatives
        self._uy = 0.0
        self._uz = 0.0
        
        self._vx = 0.0
        self._vy = 0.0
        self._vz = 0.0
        
        self._wx = 0.0
        self._wy = 0.0
        self._wz = 0.0
        
        ######### subsets ###########
        self.subset1 = None
        self.subset2 = None
        self.psubset2 = None
        
        ######### set connectivity  ###########
        self.connectList = []
        
        ######### correlation value ############
        self.correlationValue = None  
        
        ######### condition (convergence)value #############
        self.condition = False
        
        ######### for testing the optimization  ############
        self.opt_x = 0.0
        self.opt_y = 0.0
        self.opt_z = 0.0
            
    def addInitial(self, warp):
        """
           warp is [3 by 4] matrix(array)
           
            P =[[p1, p4, p7, p10]
                [p2, p5, p8, p11]
                [p3, p6, p9, p12]];
                 
                 [ux uy uz u
                  vx vy vz v
                  wx wy wz w]
    
            Convert affine warp parameters into 3 x 3 warp matrix
            NB affine parameterised as [1 + p1,    p4,     p7,     p10; [x]
                                         p2,     1 + p5,   p8,     p11; [y]
                                         p3,       p6,     p9+1,   p12  [z]
                                         0         0        0       1 ] [1]
        """
        if warp.shape[0] != 3 and warp.shape[1] != 4:
            print "warp matrix(array) input wrong size!"
        
        if not self.initialGuess:
            
            self._u = warp[0][3]
            self._v = warp[1][3]
            self._w = warp[2][3]
            
            self._ux = warp[0][0] # partial derivatives
            self._uy = warp[0][1]
            self._uz = warp[0][2]
            
            self._vx = warp[1][0]
            self._vy = warp[1][1]
            self._vz = warp[1][2]
            
            self._wx = warp[2][0]
            self._wy = warp[2][1]
            self._wz = warp[2][2]
            
            self.initialGuess = True # set the inital guess value to true
    def getNumber(self):
        return self.number
            
    def addVector(self, Z, X, Y):
        """
        this method adds a initial search vectors
         
        """     
        self._initZ = Z
        self._initX = X
        self._initY = Y
        
    def getOldCoordinator(self):
        oldz = self._z 
        oldx = self._x 
        oldy = self._y 
        return oldz, oldx, oldy
        
    def getNewCoordinator(self):
        """
        return the new coordinator after applied correlation search
        """
        newz = self._z + self._initZ 
        newx = self._x + self._initX
        newy = self._y + self._initY  
        #print self._initX

        #print newz, newx, newy
        return newz, newx, newy
    
    def getWarp(self):
        warp = np.zeros([3,4])
        warp[0][3] = self._u  
        warp[1][3] = self._v 
        warp[2][3] = self._w  
            
        warp[0][0] = self._ux   # partial derivatives
        warp[0][1] = self._uy  
        warp[0][2] = self._uz  
            
        warp[1][0] = self._vx  
        warp[1][1] = self._vy  
        warp[1][2] = self._vz  
            
        warp[2][0] = self._wx  
        warp[2][1] = self._wy  
        warp[2][2] = self._wz  
        return warp
        
    def addConnectivity(self, Node):
        self.connectList.extend(Node)
        
    def getConnectivity(self):
        return self.connectList
        
    def addCorrelationValue(self, value):
        self.correlationValue = value
        
    def getCorrelationValue(self):
        return self.correlationValue
        
    def boundarypadding(self,subset,npad):
        # pad dimension is according to (0,1,2)
        x = np.pad(subset,pad_width=npad, mode='constant', constant_values=0)
        return x
        
    def addSubset(self, sub1, sub2, psub2):
        self.subset1 = sub1
        self.subset2 = sub2
        self.psubset2 = psub2
        
    def addCorrelation(self,value):
        self.correlationValue = value
        
    def getKLTsubset(self):
        klt1 = self.subset1 
        x, y, z = klt1.shape
        x1 = x//2
        y1 = y//2
        z1 = z//2
        
        x2, y2, z2 = self.psubset2.shape
        
        centerX = (x2-1)/2 # find the center of psubset
        centerY = (y2-1)/2
        centerZ = (z2-1)/2
        
        newX = centerX + self._initZ # initZ is corresponding to X 
        newY = centerY + self._initX # initX ---> Y
        newZ = centerZ + self._initY # initY ---> Z
                                            
        klt2 = self.psubset2[newX-x1:newX+x1+1, 
                             newY-y1:newY+y1+1, 
                             newZ-z1:newZ+z1+1]
        return klt1, klt2
        
    def getSubset(self, image1, image2):
        """
        this method return two subsets regarding this node.
        """        
        coordinate1 = np.array([self._z, self._x, self._y]) 
        coordinate2 = np.array([self._z, self._x, self._y])
        subsearch = self.subSearch
        # put a boundaries in image1
        padimage1 = self.boundarypadding(image1,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
        padimage2 = self.boundarypadding(image2,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
        
        subset1 = padimage1[coordinate1[0]-self.subset:coordinate1[0]+self.subset+1,
                            coordinate1[1]-self.subset:coordinate1[1]+self.subset+1,
                            coordinate1[2]-self.subset:coordinate1[2]+self.subset+1,] # index is from 0 to n-1
                         
        subset2 = padimage2[coordinate2[0]-self.subset:coordinate2[0]+self.subset+1,
                            coordinate2[1]-self.subset:coordinate2[1]+self.subset+1,
                            coordinate2[2]-self.subset:coordinate2[2]+self.subset+1,]  

        return subset1, subset2
        
    def calculateCorrelation(self, subset1, subset2):
        result = CorrelationFunction.DCC(subset1, subset2)
        print "moving vector is: "+str(result[0])
        print "correlation is : " + str(result[1])
        self._initZ = result[0][0] # displacement in X direction
        self._initX = result[0][1]
        self._initY = result[0][2]
        return result[1] # return the minimum value
        
if __name__ == "__main__":
    
    print "in class testing!"
    
    # need to build a test unit 
    print "This is a test" 
    #import CorriFuc 
    import time
    import scipy.io as sio
    #from LoadDicom import ReadDicom
    #from DataVisu import DataVisulization
    #from PoissonDisk import DiskSamples
        
    import InputOutput 
    import SpatialSampling
    import Visulization
    import CorrelationFunction
    import scipy.ndimage as ndimage
    
    from mayavi import mlab
    import mayavi
    from scipy import spatial
    import numpy as np
    from scipy import spatial
    from scipy import linalg  
          
    #mat_contents1 = sio.loadmat('/Users/junchaowei/Desktop/Python_DVC2/UP_Research/793-crop-LC-smoothed.mat')
    #mat_contents2 = sio.loadmat('/Users/junchaowei/Desktop/Python_DVC2/UP_Research/808-crop-LC-smoothed.mat')
    
    #image1 = mat_contents1['volume']
    #image2 = mat_contents2['volume']
    
    PathDicom1 = "/Users/junchaowei/Desktop/Python_DVC2/Flex/Dicom_Flex_04_14_2016"
    PathDicom2 = "/Users/junchaowei/Desktop/Python_DVC2/Flex/Dicom_Flex_04_14_2016_2" 
    
    #PathDicom1 = "/Users/junchaowei/Desktop/Python_DVC2/CT_Data/Ball1_DICOm_ref"
    #PathDicom2 = "/Users/junchaowei/Desktop/Python_DVC2/CT_Data/Ball1_DICOm_def"
 
    dicom_volume1 = InputOutput.ReadDicom(PathDicom1)
    dicom_volume2 = InputOutput.ReadDicom(PathDicom2)
            
    dicom_volume1.loadfiles()
    dicom_volume2.loadfiles()
        
    #self.image1 = dicom_volume1.ArrayTiff
    image1 = ndimage.filters.gaussian_filter(dicom_volume1.DicArray,3)
    image2 = ndimage.filters.gaussian_filter(dicom_volume2.DicArray,3)	

    threshold1 = 1200
    threshold2 = 1200
           
    def mesh_dist(image, threshold):
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
        Samples2 = SpatialSampling.DiskSamples(image1, dist_img, threshold)
        Samples2.GenSamples(5, 80000)
        point_arr2 = Samples2.samples
        sample_point2 = np.transpose(point_arr2)
        Visulization.DataVisulization(image1, threshold).scatterplot(point_arr2)
        return sample_point2
        
    dist_img1 = mesh_dist(image1, threshold1)
    point_arr1 = sampling_def(image1, dist_img1, threshold1)
    
    dist_img2 = mesh_dist(image2, threshold2)
    point_arr2 = sampling_def(image2, dist_img2, threshold2)
    
    #sio.savemat('/Users/junchaowei/Desktop/OpenGL_Connection/PCL_point_arr1.mat', {'p':point_arr1})
    #sio.savemat('/Users/junchaowei/Desktop/OpenGL_Connection/PCL_point_arr2.mat', {'p':point_arr2})
    
    stop
    
    harrisRadius = 7
    harrisCorner = CorrelationFunction.CornerDetector(image1, harrisRadius) 
    HarrisPoint = harrisCorner.harrisPoint(point_arr1, 100)
    
    sio.savemat('/Users/junchaowei/Desktop/OpenGL_Connection/harris_nodes.mat', {'p':HarrisPoint})
    
    Visulization.DataVisulization(image1, 6000).scatterplot(HarrisPoint.T)
    
    #harrisCorner2 = CorrelationFunction.CornerDetector(image1, harrisRadius) 
    #tem = np.zeros((1,3))
    #tem[0,:] = np.array((312.0, 64.0, 111.0))
    #harrisCorner2.harrisPoint(tem, 10*10**12)
    #harrisSub = harrisCorner.visuHarris(HarrisPoint[1,:])
    
    vl = Visulization.DataVisulization(image1, 1200)
    vl.contour3d()
    
    stop0
    
    plt.figure(3)
    imshow(dist_img1[:,:,50])
    show()
    plt.figure(5)
    imshow(image1[:,:,50])
    show()
    
    sio.savemat('/Users/junchaowei/Desktop/OpenGL_Connection/surf_nodes.mat', {'p':vl.surf_points})
    sio.savemat('/Users/junchaowei/Desktop/OpenGL_Connection/surf_elements.mat', {'p':vl.surf_elements})


    harris_point_arr1 = np.array(HarrisPoint)
    
    mesh2 = SpatialSampling.MeshDelaunay(point_arr1.T) # use the original point_arr1 
    mesh2.mesh()
    mesh2.alpha_shape(25)
    mesh2.mesh_triangle() # find unique triangles
    
    for dummy_i in range(1):   
        Mesh_optimizer = SpatialSampling.MeshOptimizer(point_arr1.T,mesh2.ntri)
        Mesh_optimizer.edge_connector_smooth()
        
    #mesh2.points
    
    ###
    #image2_3, point_arr2_3 = CorrelationFunction.ICP(image2, point_arr1, point_arr2,5)
    ###
    stop
    
    vl = Visulization.DataVisulization(image1, 1200)
    vl.contour3d()
    vl = Visulization.DataVisulization(image2, 1200)
    vl.contour3d()    
                
    NodeList = []
    subSize = 25
    subSearch = 13  
    dictions = CorrelationFunction.calculateSubset(image1, image2, mesh2.points, subSize, subSearch)
    sub1 = dictions[0] 
    nodesDiction={}

    for i in range(mesh2.points.shape[0]):
        print i
        # initializing a instance 
        tempNode = Node([mesh2.points[i][0], mesh2.points[i][1], mesh2.points[i][2]], subSize, subSearch, i)
        
        tempNode.addSubset(dictions[0].get(i), dictions[1].get(i), dictions[4].get(i)) 
        # the output index dictions[index]: 0-subset1, 1-subset2, 2-correlation value, 3-displacements, 4-psubset2
        # add the correlation value to each node
        tempNode.addCorrelation(dictions[2].get(i))
        # add initial search vector to each node
        tempNode.addVector(dictions[3].get(i)[0], dictions[3].get(i)[1], dictions[3].get(i)[2]) # from dictionary adding three displacements (z,x,y) value into a node instance
        
        # calculating connecting edge point to a point
        connected_points = CorrelationFunction.ConnectMesh(mesh2.ntri,i)
        
        # add the connectivity to the node
        tempNode.addConnectivity(connected_points)
        
        # append the node to the list
        NodeList.append(tempNode)
        
        # put node in the dictionary for later search
        nodesDiction[i] = tempNode
        
    # sort the node according to the correlation value
    NodeList.sort(key=lambda x: x.correlationValue, reverse=False) # sort the list in place
    
    index = 0
    outputPoints = []
    outputColor = []
    size = str(len(NodeList))
    
    for eachNote in NodeList:
        
        print "the computation process is: "+str(NodeList.index(eachNote))+ " / "+ size
        print "the correlation value of node is: "+ str(eachNote.correlationValue)
        index += 1
        if eachNote.initialGuess == False: # if the node has not been calculated.
            
            #eachNote.subset1
            #eachNote.subset2
    
            subsets = eachNote.getKLTsubset()
    
            #print NodeList[nodei].getOldCoordinator()
            #print NodeList[nodei].getNewCoordinator()
    
            old_loc = eachNote.getOldCoordinator()
            new_loc = eachNote.getNewCoordinator()
            
            image1_coordinator = np.ceil(np.array([old_loc[0], old_loc[1], old_loc[2]]))
            image2_coordinator = np.ceil(np.array([new_loc[0], new_loc[1], new_loc[2]]))   
            
            connected = eachNote.getConnectivity()
            
            iniWarp = np.zeros([3,4])
            #count = 0.0
            #for item in connected:
            #    if nodesDiction.get(item).initialGuess != False:
            #        count +=1 # counting how many nodes connected.
            #        iniWarp = iniWarp + nodesDiction.get(item).getWarp() # sum all the warp value into a matrix
                    
            #if count != 0: # if there is a connection for the node:
            #    iniWarp = iniWarp/count
            
            try:
                good = CorrelationFunction.KLT(image1, image2, image1_coordinator, image2_coordinator, subSize, iniWarp)
            except ValueError:
                print("jump the loop: "+str(index))
                continue
        
            translation = good[1]
            if np.sqrt(sum(translation**2))<5:
                eachNote.addInitial(good[0])
                eachNote._initZ = eachNote._initZ + translation[0]
                eachNote._initX = eachNote._initX + translation[1]
                eachNote._initY = eachNote._initY + translation[2]
                
                eachNote.opt_x = translation[0]
                eachNote.opt_y = translation[1]
                eachNote.opt_z = translation[2]
                
                eachNote.condition = good[4]
            else:
                print "large translation values!!!"
            
            outputPoints.append(eachNote.getOldCoordinator())
            outputColor.append(eachNote.condition)
    
    outputPoints = np.array(outputPoints)
#    sio.savemat('/Users/junchaowei/Desktop/OpenGL_Connection/harris_nodes_reshape.mat', {'p':outputPoints})
#    sio.savemat('/Users/junchaowei/Desktop/OpenGL_Connection/harris_nodes_reshape_color.mat', {'c':outputColor})
            
    stop2
    
    subsets = eachNote.getKLTsubset()
        
#    if index == 9:
#       break
        
    subset1 = subsets[0]
    subset2 = subsets[1]
    
    vl = Visulization.DataVisulization(subset1, 1200)
    vl.contour3d()
    vl = Visulization.DataVisulization(subset2, 1300)
    vl.contour3d()
    
    vl = Visulization.DataVisulization(eachNote.psubset2, 110)
    vl.contour3d()
        
    #valueCorrelation = np.zeros(mesh2.points.shape[0])
            
    #sio.savemat('/Users/junchaowei/Desktop/sub1.mat', {'v':sub1}) 
    #sio.savemat('/Users/junchaowei/Desktop/sub2.mat', {'v':sub2}) 
    #sio.savemat('/Users/junchaowei/Desktop/psub2.mat', {'v':NodeList[nodei].psubset2})
    
    t1 = good[2]
    t2 = good[3]

    vl = Visulization.DataVisulization(ndimage.sobel(image1,0), ndimage.sobel(image1,0).mean())
    vl.contour3d()
    vl = Visulization.DataVisulization(image2, 1400)
    vl.contour3d()
    
    dispx = np.zeros(len(nodesDiction))
    dispy = np.zeros(len(nodesDiction))
    dispz = np.zeros(len(nodesDiction))
    nodes_dis = []
    for icolorNode in nodesDiction.keys():
        dispx[icolorNode] = nodesDiction.get(icolorNode).opt_z
        dispy[icolorNode] = nodesDiction.get(icolorNode).opt_x
        dispz[icolorNode] = nodesDiction.get(icolorNode).opt_y
        nodes_dis.append(nodesDiction.get(icolorNode).getOldCoordinator())
        
    colors = dispz
    @mlab.show
    def main():
        mesh2.view(mesh2.polydata(mesh2.points.T,colors))  
    main()
    mayavi.mlab.colorbar(object=None, title="U(Z)(voxel)", orientation="vertical", nb_labels=20, nb_colors=None, label_fmt=None)
    "colorbar reference: http://docs.enthought.com/mayavi/mayavi/auto/mlab_decorations.html"
    
    from scipy import interpolate
    inter_z = interpolate.LinearNDInterpolator(nodes_dis, dispz, fill_value=0.0, rescale=True)
    print "imported!3"
   
   
    plot = Visulization.DataVisulization(image1,1300)
    plot.contour3d()
    surf = plot.surf_points
    
    extrapo_z = inter_z.__call__(surf)
    plot.plot3dpoint(extrapo_z)
    
    #plot = Visulization.DataVisulization(image1,image1.mean())
    #plot.contour3d()
    
    #########################  Strains ##############################
    
    ################  smooth the neigbour displacement point using the neighbor average #########
    edge1 = np.sort(np.array([mesh2.ntri[:,0],mesh2.ntri[:,1]]),axis = 0)
    edge2 = np.sort(np.array([mesh2.ntri[:,0],mesh2.ntri[:,2]]),axis = 0)
    edge3 = np.sort(np.array([mesh2.ntri[:,0],mesh2.ntri[:,3]]),axis = 0)
    edge4 = np.sort(np.array([mesh2.ntri[:,1],mesh2.ntri[:,2]]),axis = 0)
    edge5 = np.sort(np.array([mesh2.ntri[:,1],mesh2.ntri[:,3]]),axis = 0)
    edge6 = np.sort(np.array([mesh2.ntri[:,2],mesh2.ntri[:,3]]),axis = 0)
    edge = np.concatenate((edge1,edge2,edge3,edge4,edge5,edge6),axis=1)
    trans_edge = edge.T
    b = np.ascontiguousarray(trans_edge).view(np.dtype((np.void, trans_edge.dtype.itemsize * trans_edge.shape[1])))
    _, idx = np.unique(b, return_index=True)
    mesh2.edge = trans_edge[idx] # this is the unique edges

    target_points = mesh2.points  # target_points shold be fixed at each time
    
    xdis = dispx
    ydis = dispy
    zdis = dispz
    
    gradi_xdis = np.zeros_like(xdis)
    gradi_ydis = np.zeros_like(ydis)   
    gradi_zdis = np.zeros_like(zdis)                                                                                                                                                     
    for i in range(mesh2.points.shape[0]):
        print i
        connector1 = mesh2.edge[mesh2.edge[:,0] == i,1] # find the connectivity from the first column of edge
        connector2 = mesh2.edge[mesh2.edge[:,1] == i,0] # find the connectivity from the second column of edge
        
        if connector1.size !=0  and connector2.size != 0:
            connector = np.concatenate((connector1,connector2))
        elif connector1.size == 0:
            connector = connector2
        elif connector2.size == 0:
            connector = connector1
        #connector = np.concatenate((connector1,connector2),axis=1)
        
        connecting_points = target_points[connector,:]
        connecting_xdis = xdis[connector]
        connecting_ydis = ydis[connector]
        connecting_zdis = zdis[connector]
    
        if connector.shape[0] != 0:
        ##############  uniform laplacian #############
            #print connector
            vector = np.sum(connecting_points - target_points[i,:], axis = 0)/connector.shape[0] 
            direction = np.sign(connecting_points[:,0] - target_points[i,0])
            
            # for small area the length of bar in each direction can be approximated by the length of bar
            gradi_norm = np.sqrt(np.sum((connecting_points - target_points[i,:])**2, axis = 1))
            
            #gradi_norm2 = np.zeros((len(gradi_norm),1))
            #gradi_norm2[:,0]=gradi_norm
            #gradi_vect = (connecting_points - target_points[i,:])/gradi_norm2
            
            #lenth_gradi_vect1 = abs(connecting_points[:,0] - target_points[i,0])
            #lenth_gradi_vect2 = abs(gradi_vect[:,1]*gradi_norm)
            #lenth_gradi_vect3 = abs(gradi_vect[:,2]*gradi_norm)
            direction_diff = (connecting_zdis - ydis[i])*direction
            # the transpose .T[0] keep the same shape of (connecting_xdis - xdis[i]) with the gradi_norm 
            #gradi_xdis[i] = np.sum(direction_diff.T[0]/gradi_norm, axis = 0)/connector.shape[0] 
            #gradi_ydis[i] = np.sum(direction_diff.T[0]/gradi_norm, axis = 0)/connector.shape[0] 
            gradi_zdis[i] = np.sum(direction_diff.T[0]/gradi_norm, axis = 0)/connector.shape[0] 
            
            #gradi_zdis[i] = np.sum(direction_diff.T[0]/gradi_norm, axis = 0)/connector.shape[0] 
            #gradi_ydis[i] = np.sum((connecting_ydis - extrapo_y[i]).T[0]/gradi_norm, axis = 0)/connector.shape[0] 
            #gradi_zdis[i] = np.sum((connecting_zdis - extrapo_z[i]).T[0]/gradi_norm, axis = 0)/connector.shape[0] 
    
    colors = gradi_zdis
    
    @mlab.show
    def main():
        mesh2.view(mesh2.polydata(Mesh_optimizer.points.T,colors))  
    main()
    mayavi.mlab.colorbar(object=None, title="Strain(Z)(voxel)", orientation="vertical", nb_labels=10, nb_colors=None, label_fmt=None)
    "colorbar reference: http://docs.enthought.com/mayavi/mayavi/auto/mlab_decorations.html"
            
            