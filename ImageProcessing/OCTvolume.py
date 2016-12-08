"""
this creates volume object

"""

import os
import ImageProcessing
import SpatialSampling
import CorrelationFunction
import numpy as np
import scipy.io as sio

class OCTimage(object):
    """
    root class OCT image
    """ 
    def __init__(self, info, start=None, end=None):
        
        self._start = start
        self._end = end
        
        if os.path.isdir(info):
            self.address = info # address to access the tiff images
            self._volume_data = self.readTiff() 
               
        elif os.path.isfile(info):
            print "read .mat file!"
            self.address = info
            self._volume_data = self.readMat()
        else:
            raise TypeError("the input address must be a string file name or directory!")
     
        if start or end:
            self._state = True # False if none cropped (default), True Ohterwise
        else:
            self._state = False
        
    def readTiff(self):
        if self._start and self._end:
            out = ImageProcessing.readPartialOCTstack(self.address, self._start, self._end)
            print "read partial the image stacks: " + str(self._start)+ " to " + str(self._end)
            print "read the image stack size: " + str(out.shape)     
    
        else:
            out = ImageProcessing.readOCTstack(self.address)
            print "read all the image stacks"
            print "read the image stack size: " + str(out.shape)
            
        ##############  output ###############
        print "the volume max: " + str(out.max())
        print "the volume min: " + str(out.min())
        print "the volume mean: " + str(out.mean())    
        return out
        
    def readMat(self):
        volume000 = sio.loadmat(self.address) 
        for i in range(len(volume000)):
            if type(volume000.get(volume000.keys()[i])) == np.ndarray:
                return volume000.get(volume000.keys()[i])                                                                                                                                               

class RefVolume(OCTimage):
    """
    reference volume class:
        this class should read the tiff image stack and convert it into nd-arrary
        
        it should contain several class methods including:
        1) create volume from images
        2) noise reduction volume, dynamic range, ostu method
        3) create binary domain
        4) volume discretization 
        ...        
    """
    def __init__(self, info,start=None,end=None):
        """
        initialize the class
        """
        OCTimage.__init__(self, info, start, end)
        self._convex_xy = None # the list stores the x dimension of convex, the z dimension is as the index of element in list
        self._periPoint = None
        self._logicSlice = None
        self._block_volume = None
        self._denoise_volume = None
        self._binary_volume = None

    def denoise2D(self,volume, size, inRange):
        """
        deprecated@@!
        median_filter on the block volume, if size = 0, the median filter is ignored
        dynamics range removes the background salt and pepper noise
        size --- for median filter
        inRange --- for Dynamic Range 
        """
        print "using 2D median filter!"
        if size == 0:
            print "no median filter!"
            self._denoise_volume = ImageProcessing.DynamicRangeImage(volume, inRange)
        else: 
            temp_volume = ImageProcessing.medOCT(volume, size)
            #temp_volume = ndimage.gaussian_filter(temp_volume,1)
            self._denoise_volume = ImageProcessing.DynamicRangeImage(temp_volume, inRange)
    
    def denoise3D(self, volume, size, inRange):
        """
        this uses 3D median filter
        median_filter on the block volume, if size = 0, the median filter is ignored
        dynamics range removes the background salt and pepper noise
        size --- for median filter
        inRange --- for Dynamic Range 
        """
        print "using 3D median filter!"
        if size ==0:
            print "no median filter!"
            self._denoise_volume = ImageProcessing.DynamicRangeImage(volume, inRange)
        else:
            temp_volume = ImageProcessing.SpaceMedianFilter(volume, size)
            #temp_volume = ndimage.gaussian_filter(temp_volume,1)
            self._denoise_volume = ImageProcessing.DynamicRangeImage(temp_volume, inRange)
    
    def simpleBinary(self, block_size):
        self._binary_volume = ImageProcessing.sliceThreshold(self._block_volume, block_size=block_size)
  
    def points(self, number = None):
        """
        selection convex points in the image
        number --- the slice number in the volume
        """
        if number:
            PointBuilder = ImageProcessing.KeyPoints(self._volume_data[:,:,number])
        else:
            sx, sy, sz = self._volume_data.shape
            PointBuilder = ImageProcessing.KeyPoints(self._volume_data[:,:,sz//2])
                    
        return PointBuilder
        
    def segmentation(self, points, useSnake = True): 
        
        ## make ellipse
        if len(points)==0:
            print "the input points size is ZERO!"
            return None
        
        Initial_ellipsePoints = ImageProcessing.FindEllipse(points)
        
        self._convex_xy = Initial_ellipsePoints # back up the data
                
        outPoints = []
        
        if self._state == False:
            for i in range(self._volume_data.shape[2]): # iterating through the vertical stacks
                print i
                slice_snake = np.empty([Initial_ellipsePoints.shape[0],3]) # Make N by 3 points
                
                if useSnake == True:
                    snake = ImageProcessing.ActiveEllipse(self._volume_data[:,:,i], Initial_ellipsePoints)
                else:
                    snake = Initial_ellipsePoints
                    
                slice_snake[:,0] = snake[:,0] # create 3d points
                slice_snake[:,1] = snake[:,1]
                slice_snake[:,2] = i
                outPoints.append(slice_snake)
                
        elif self._state == True:
            for i in range(self._start, self._end+1): # iterating through the vertical stacks
                print i
                slice_snake = np.empty([Initial_ellipsePoints.shape[0],3]) # Make N by 3 points 
                
                if useSnake == True:
                    snake = ImageProcessing.ActiveEllipse(self._volume_data[:,:,i], Initial_ellipsePoints)
                else:
                    snake = Initial_ellipsePoints
                    
                slice_snake[:,0] = snake[:,0] # create 3d points
                slice_snake[:,1] = snake[:,1]
                slice_snake[:,2] = i
                outPoints.append(slice_snake)
                
        for i in range(len(outPoints)): # reshape the list
            if i == 0: fb = outPoints[i]  
            else: fb = np.concatenate((fb,outPoints[i]),axis=0)
            
        self._periPoint = fb
        
    def blockImage(self,volume,points):
        """
        this use the perPoint to create convex hull and logical operation for each slice
        """
        self._logicSlice = ImageProcessing.LogicRegion(self._denoise_volume[:,:,self._denoise_volume.shape[2]//2], self._convex_xy)
        self._block_volume = np.empty_like(self._denoise_volume)
        if self._state == False:
            for i in range(self._volume_data.shape[2]):
                self._block_volume[:,:,i] = self._denoise_volume[:,:,i] * self._logicSlice
                
        elif self._state == True:
             for i in range(self._start, self._end+1): 
                self._block_volume[:,:,i] = self._denoise_volume[:,:,i] * self._logicSlice
        
    def threshold(self, block_size):
        return ImageProcessing.sliceThreshold(self._volume_data, block_size)
        
    def getOriginVolume(self):
        return self._volume_data  
        
    def getBlockVolume(self):
        return self._block_volume   
        
    def getBinaryVolume(self):
        return self._binary_volume
        
    def getDenoiseVolume(self):
        return self._denoise_volume
                
class DefVolume(OCTimage):
    """
    reference volume class:
        this class should read the tiff image stack and convert it into nd-arrary
        
        it should contain several class methods including:
        1) create volume from images
        2) noise reduction volume, dynamic range, ostu method
        3) create binary domain
        4) volume discretization 
        ...        
    """
    def __init__(self, address,start=None,end=None):
        """
        initialize the class, deformed volume do not have the mesh sampling
        
        """
        OCTimage.__init__(self, address, start, end)
        self._convex_xy = None # the list stores the x dimension of convex, the z dimension is as the index of element in list
        self._periPoint = None
        self._logicSlice = None
        self._block_volume = None # logic selected region
        self._denoise_volume = None
        self._binary_volume = None
        
    def denoise2D(self,volume, size,inRange):
        """
        deprecated@@!
        this uses 2D median filter
        median_filter on the block volume, if size = 0, the median filter is ignored
        dynamics range removes the background salt and pepper noise
        size --- for median filter
        inRange --- for Dynamic Range 
        """
        print "using 2D mediam filter!"
        if size == 0:
            print "no median filter!"
            self._denoise_volume = ImageProcessing.DynamicRangeImage(volume, inRange)
        else: 
            temp_volume = ImageProcessing.medOCT(volume, size)
            #temp_volume = ndimage.gaussian_filter(temp_volume,3)
            self._denoise_volume = ImageProcessing.DynamicRangeImage(temp_volume, inRange)
            
    def denoise3D(self,volume, size, inRange):
        """
        this uses 3D median filter
        median_filter on the block volume, if size = 0, the median filter is ignored
        dynamics range removes the background salt and pepper noise
        size --- for median filter
        inRange --- for Dynamic Range 
        """
        print "using 3D median filter!"
        if size ==0:
            print "no median filter!"
            self._denoise_volume = ImageProcessing.DynamicRangeImage(volume, inRange)
        else:
            temp_volume = ImageProcessing.SpaceMedianFilter(volume, size)
            #temp_volume = ndimage.gaussian_filter(temp_volume,3)
            self._denoise_volume = ImageProcessing.DynamicRangeImage(temp_volume, inRange)
        
    def simpleBinary(self, volume, block_size):
        self._binary_volume = ImageProcessing.sliceThreshold(volume, block_size=block_size)
        
    def points(self, number = None):
        """
        selection convex points in the image
        number --- the slice number in the volume
        """
        if number:
            PointBuilder = ImageProcessing.KeyPoints(self._volume_data[:,:,number])
        else:
            sx, sy, sz = self._volume_data.shape
            PointBuilder = ImageProcessing.KeyPoints(self._volume_data[:,:,sz//2])
                    
        return PointBuilder
        
    def segmentation(self, points, useSnake = True): 
        
        ## make ellipse
        if len(points)==0:
            print "the input points size is ZERO!"
            return None
        
        Initial_ellipsePoints = ImageProcessing.FindEllipse(points)
        
        self._convex_xy = Initial_ellipsePoints # back up the data

        outPoints = []
        
        if self._state == False:
            for i in range(self._volume_data.shape[2]): # iterating through the vertical stacks
                print "uncropped: " + str(i)
                slice_snake = np.empty([Initial_ellipsePoints.shape[0],3]) # Make N by 3 points
                
                if useSnake == True:
                    snake = ImageProcessing.ActiveEllipse(self._volume_data[:,:,i], Initial_ellipsePoints)
                else:
                    snake = Initial_ellipsePoints
                    
                slice_snake[:,0] = snake[:,0] # create 3d points
                slice_snake[:,1] = snake[:,1]
                slice_snake[:,2] = i
                outPoints.append(slice_snake)
                
        elif self._state == True:
            for i in range(self._start, self._end+1): # iterating through the vertical stacks
                print "cropped: " + str(i)
                slice_snake = np.empty([Initial_ellipsePoints.shape[0],3]) # Make N by 3 points 
                
                if useSnake == True:
                    snake = ImageProcessing.ActiveEllipse(self._volume_data[:,:,i], Initial_ellipsePoints)
                else:
                    snake = Initial_ellipsePoints
                    
                slice_snake[:,0] = snake[:,0] # create 3d points
                slice_snake[:,1] = snake[:,1]
                slice_snake[:,2] = i
                outPoints.append(slice_snake)
                
        for i in range(len(outPoints)): # reshape the list
            if i == 0: fb = outPoints[i]  
            else: fb = np.concatenate((fb,outPoints[i]),axis=0)
            
        self._periPoint = fb
        
    def blockImage(self):
        """
        this use the perPoint to create convex hull and logical operation for each slice
        """
        self._logicSlice = ImageProcessing.LogicRegion(self._denoise_volume[:,:,self._denoise_volume.shape[2]//2], self._convex_xy)
        self._block_volume = np.empty_like(self._denoise_volume)
        if self._state == False:
            for i in range(self._denoise_volume.shape[2]):
                self._block_volume[:,:,i] = self._denoise_volume[:,:,i] * self._logicSlice
                
        elif self._state == True:
             for i in range(self._start, self._end+1): 
                self._block_volume[:,:,i] = self._denoise_volume[:,:,i] * self._logicSlice

    def threshold(self, block_size):
        return ImageProcessing.sliceThreshold(self._volume_data, block_size)
           
    def getOriginVolume(self):
        return self._volume_data  
        
    def getBlockVolume(self):
        return self._block_volume  
        
    def getBinaryVolume(self):
        return self._binary_volume
        
    def getDenoiseVolume(self):
        return self._denoise_volume
    
if __name__ == "__main__":
    print "this is a unit test!"
    import numpy as np
    import matplotlib.pyplot as plt # volume.getVolume().shape[2]//2
    from skimage.segmentation import active_contour
    import Visualization
    from scipy import ndimage
    
    #Refvolume2 = RefVolume("/Users/junchaowei/Desktop/Python_DVC2/UP_Research/WholeRegionRealData/JunchaoFirstRun/M6_OD/M6_OD_baseline_C-scan/")
    #Defvolume = DefVolume("/Users/junchaowei/Desktop/Python_DVC2/UP_Research/WholeRegionRealData/JunchaoFirstRun/M6_OD/M6_OD_125_C-scan/")

    Refvolume = RefVolume("/Users/junchaowei/Desktop/SpaceRegistration_000_125/volume000.mat")
    Defvolume = DefVolume("/Users/junchaowei/Desktop/SpaceRegistration_000_125/volume125_regi.mat") 
    
   # ra = Refvolume.points()
   # da = Defvolume.points()
   #
   # rb = Refvolume.segmentation(np.array([ra.xs[1:],ra.ys[1:]]).T, useSnake=True)
   # db = Defvolume.segmentation(np.array([da.xs[1:],da.ys[1:]]).T, useSnake=True)

    Refvolume.denoise2D(Refvolume.getOriginVolume(), 5, (50,255))
    Defvolume.denoise2D(Defvolume.getOriginVolume(), 5, (50,255)) 

    #Refvolume.simpleBinary(101)
    #Defvolume.simpleBinary(101)
   
    de_volume1 = Refvolume.getDenoiseVolume()
    de_volume2 = Defvolume.getDenoiseVolume()

    #plt.figure(4)
    #plt.imshow(Refvolume2._volume_data[:,:,-10],interpolation='nearest',cmap='Greys_r')
    #plt.show()
    
#            
#    points_volume1 = Refvolume.getBinaryVolume()
#    points_volume2 = Defvolume.getBinaryVolume() 
#            
#    level1 = 0.5
#    level2 = 0.5
#            
#    number1 = 100000 # base line
#    number2 = 100000
#        
#    sampling = SpatialSampling.RandomSampling
#    icp = CorrelationFunction.ICP
#        
#        #modi_volume = orig_volume2
#        #for i in range(iterations):
#        #print "iteration #: " + str(i)
#        
#    volume_points_1 = sampling(ndimage.gaussian_filter(points_volume1,3), level1, number1) # generating new point avoiding local minimum
#    volume_points_2 = sampling(ndimage.gaussian_filter(points_volume2,3), level2, number2)
#
#    new_volume2, change_points, rotation_matrix_list, translation_vector_list = icp(orig_volume2,volume_points_1,volume_points_2,5)
#            
    plot = Visualization.DataVisulization(de_volume1, 0.3)
    plot.contour3d()
    plot = Visualization.DataVisulization(de_volume2, 0.3)
    plot.contour3d()
    
    test = ndimage.median_filter(Refvolume.getOriginVolume()[:,:,60], 10)
    
    plt.figure(1)
    plt.imshow(test,cmap='Greys_r')
    plt.show()
    
    from skimage import exposure
    
    denoise_volume = exposure.rescale_intensity(test, in_range=(50,255))
    
    plt.figure(2)
    plt.imshow(denoise_volume,cmap='Greys_r')
    plt.show()
    
    denoise_volume2 = exposure.equalize_adapthist(denoise_volume,clip_limit=0.1)
    
    plt.figure(3)         
    plt.imshow(ndimage.gaussian_filter(denoise_volume2,2),cmap='Greys_r')
    plt.show()   
    
    from skimage.filters import threshold_otsu, threshold_adaptive
    denoise_volume3 = threshold_adaptive(ndimage.gaussian_filter(denoise_volume2,2), 5, offset=0)
    
    plt.figure(4)
    plt.imshow(denoise_volume3,cmap='Greys_r')
    plt.show()
    
#        
#    plot = Visualization.DataVisulization(ndimage.gaussian_filter(new_volume2,2), 80)
#    plot.contour3d()
#   
#    stop
#    

    
#    #plt.imshow(active[:,:,50],cmap='Greys_r')
#    #plt.show()
#    plt.imshow(ndimage.gaussian_filter(Refvolume._binary_volume,2)[:,:,50],cmap='Greys_r')
#    plt.show()
#    
#    plot = Visualization.DataVisulization(ndimage.gaussian_filter(Refvolume._binary_volume,2), 0.3)
#    plot.contour3d()
#    
#    plot = Visualization.DataVisulization(ndimage.gaussian_filter(Refvolume._block_volume,2), 80)
#    plot.contour3d()
#        
#    input_volume = ndimage.gaussian_filter(Refvolume._binary_volume,3)
#    
#    new = np.concatenate((volume_points_1.T,volume_points_2.T),axis=1)
#    Visualization.DataVisulization(Refvolume, 1).scatterplot(new)
#    Visualization.DataVisulization(Refvolume, 1).scatterplot(volume_points_2.T)
# 
    #points = np.empty([len(b)*b[0].shape[0],3]) # initlizing the point array, each slice points number is b[0].shape, the slice number is len(b)
            
    #Visualization.DataVisulization(Refvolume, 1).scatterplot(fb.T)
    
    #snake = active_contour(volume.getVolume()[:,:,40], b, w_edge=-1, w_line=1, gamma=0.001, beta=10) 
    #plt.plot(b[:,0],b[:,1])
    #plt.imshow(volume.getVolume()[:,:,40],cmap='Greys_r')
    #plt.show()
    
    #
    #threshold = volume.threshold(251)
    #import Visualization
    ##plot = Visualization.DataVisulization(ndimage.gaussian_filter(threshold,1), 0.3)
    ##plot.contour3d()
    #import matplotlib.pyplot as plt
    #plt.imshow(volume.getVolume()[:,:,30],cmap='Greys_r')
    #plt.imshow(threshold[:,:,30],cmap='Greys_r')
    #plt.show()
