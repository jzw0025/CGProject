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
    
    PathDicom1 = "Dir1"
    PathDicom2 = "Dir2" 
 
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
        
    harrisRadius = 7
    harrisCorner = CorrelationFunction.CornerDetector(image1, harrisRadius) 
    HarrisPoint = harrisCorner.harrisPoint(point_arr1, 100)
    
    Visulization.DataVisulization(image1, 6000).scatterplot(HarrisPoint.T)
    
    vl = Visulization.DataVisulization(image1, 1200)
    vl.contour3d()
  
 