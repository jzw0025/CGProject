"""
This module creates a application for the Poission Disk Sampling, and it discretizes the domain
of field. The DiskSample class creates istance for the Poisson Disk Sampler

The Ostu's Threshold Method could be applied to discriminate the different image regions.

"""
sample_key = 1
try:
    import DiskSampling as ds1
except ImportError:
    #raise ImportError('The Python-Fortran Module has not been found! Using the pure Python Module instead!')
    print 'Sampling Function, The Python-Fortran Module has not been found! Using the pure Python Module instead!'
    #import DiskSampling_PY as ds1
    import DiskSampling_PY_KDtree as ds1
    sample_key = 2
    

import numpy as np
import random
import math
from scipy.spatial import Delaunay

class DiskSamples:
    """
    This class DiskSamples creates the Poisson Disk Sampling:
       1) init_surface creates surface random points for specific ISOsurface
       2) dis_check the survey greyscale voxel, which should be within the range of distance function
       3) in_hull checks the sample to see if it is within the convex hull, which created by Delaunay triangulation
       4) GenSamples creates poisson disk samples within specified ISO value.
    """
    def __init__(self,segimag, meshsize, threshold):
        self.seg_arr = segimag
        self.samples = []
        self.sizing = meshsize
        self.threshold = threshold
        self.surf_points = [] # this array shape is (N,3)
        self.surf_elements = [] # this array shape is (N,4) ploygon
        self.hist_count = []
        self.hist_bin_edges = []
        self.level = []
        self.active_size = [0,0] # initialize the array for active sample size, which allows to use active_size[-2] and active_size[-1]
        self.mesh_scale_factor = 2 # >0, 1 =  no scale effect,  (0,1) shrinkage effect
        print "new"
        
    def getSamples(self):
        return self.samples
        
    def init_surface(self,number_points):
        """
        this method extracts the isosurface points as initial points array
        the input 'number_points' denotes the uniformly selected points from
        surface array.
        
        """
        plot = DataVisulization(self.seg_arr,self.threshold)
        plot.contour3d()
        surf = plot.surf_points
        elements = plot.surf_elements
        xsize, ysize = surf.shape
        
        # uniform randomly picking the points from the array        
        index = np.random.random_integers(0, xsize-1, number_points)
        # take the row out of numpy array by index numpy array
        self.surf_points = surf[index,:]
        self.surf_elements = elements
        
        # be careful about the shape input, the surf_point.T shape is (3,N), which is transposed from surf_point
        # visulization of random surface point
        # DataVisu.DataVisulization(self.seg_arr,self.threshold).scatterplot(self.surf_points.T)
        
    def dis_check(self,points_array,threshold):
        
        # points_array input is (3,N) numpy array
        xdim, ydim =  points_array.shape
        
        keep_points = np.array([[],[],[]])
        
        # the first point is orignal sample point, which already in the sample list, jumping over it.
        for i in xrange(1,ydim):
            
            xq = int(points_array[0,i])
            yq = int(points_array[1,i])
            zq = int(points_array[2,i])
            
            # check distance for each point in the array
            if self.sizing[xq][yq][zq] > 0 and xq!=0 and yq !=0 and zq!=0: # the distance >0, the point is inside of shape
                #print self.sizing[xq][yq][zq]
                
                new_points_array = np.array([[points_array[0,i]],[points_array[1,i]],[points_array[2,i]]])   
                keep_points = np.append(keep_points,new_points_array,axis=1)
                
        return keep_points
        
        
    def in_hull(self, p, hull):
        """
        Test if points in `p` are in `hull`
        `p` should be a (N,3) `NxK` coordinates of `N` points in `K` dimensions
        `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
        coordinates of `M` points in `K`dimensions (M,3) for which Delaunay triangulation
        will be computed
        This method is cited from overflow by user 'Juh_'
        """
        if not isinstance(hull,Delaunay):
           hull = Delaunay(hull)
           return hull.find_simplex(p) >= 0 
           # It returns a boolean array where True values indicate points that lie in the given convex hull. 
           # It can be used like this:
           # tested = np.random.rand(20,3)
           # cloud  = np.random.rand(50,3)
           # print in_hull(tested,cloud)
           
    def image_histogram(self,image,image_format):
        
        # the image format gives the bit of image, which are such as unsigned '16' or '32'
        # it needs to reshape the 3D array into a flattened array
        # using the numpy.reshape and numpy.histogram
        # numpy.rehshape : http://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
        # numpy.histogram : http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html
        
        imagex, imagey, imagez = image.shape
        number_voxels = imagex*imagey*imagez
        reshaped_image = image.reshape((1,number_voxels))
        max_int = 2**image_format
        hist, bin_edges = np.histogram(reshaped_image, bins=np.arange(max_int), density=False)
        # density =  False will give the counts of histogram while the True will give the density
        self.hist_count = hist
        self.hist_bin_edges = bin_edges
        
        # the following plots the histogram
        #import matplotlib.pyplot as plt
        #width = 0.7 * (bin_edges[1] - bin_edges[0])
        #center = (bin_edges[:-1] + bin_edges[1:]) / 2
        #plt.bar(center, hist, align='center', width=width)
        #plt.show()
        # The object-oriented interface is also straightforward:
        #fig, ax = plt.subplots()
        #ax.bar(center, hist, align='center', width=width)
        #fig.savefig("1.png")
          
    def otsu_method(self):
        """
        this method implements the Otsu's threshold method to maximize the variance between group and minimize the 
        variance within group
        The reference code is from 
        Wiki: https://en.wikipedia.org/wiki/Otsu%27s_method
        Slices: http://courses.cs.tau.ac.il/pyProg/1112a/recitations/ta10.pdf
        (This method needs to calculate the image histogram in advance)
        
        """
        hist = self.hist_count
        # this is total voxel in the image
        total = sum(self.hist_count)  
        # initialize the values
        sum_back, w_back, w_fore, var_max, threshold1, threshold2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        
        max_voxel = self.hist_bin_edges.max()
        # print max_voxel
        sum_all = 0.0
        # sum the values of all background pixels
        for t in range(max_voxel):
            sum_all += t * hist[t]
            
        for dummy_t in range(max_voxel):
            # update weights
            w_back += hist[dummy_t]
            if (w_back == 0): 
                continue
            w_fore = total - w_back
            if (w_fore == 0): 
                break
            # calculate classes means
            sum_back += dummy_t * hist[dummy_t]
            mean_back = sum_back / w_back
            mean_fore = (sum_all - sum_back) / w_fore
            # Calculate Between Class Variance
            var_between = w_back * w_fore * (mean_back - mean_fore)**2
            if (var_between >= var_max):
                threshold1 = dummy_t
                if (var_between > var_max):
                    threshold2 = dummy_t
                var_max = var_between
        self.level = (threshold1 + threshold2)/2
        
    def ostu_method_local_partition(self,image,box_size):
        """
        Create a input and output data, so it can be used externally.
        '//' is floor division operator
        '%' is the rem division operator
        e.g. 5//2 = 2 while 5%2 =1 
        """
        sx,sy,sz = image.shape
        
        pass
     
    @staticmethod # static method without object   
    def RandomSampling(image, threshold, number, meshsize=None):
        """
        random sample generation in the given domain, meshsize.
        docs: https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.uniform.html
        image --- input 3d data array
        meshsize --- input distance domain (optional)
        threshold --- surface value envoloped domain
        number --- the size/density of samples
        """
        
        sx, sy, sz = image.shape
        
        samples_x = np.random.uniform(0, sx-1, number)
        samples_y = np.random.uniform(0, sy-1, number)
        samples_z = np.random.uniform(0, sz-1, number)
        
        keep_points = []
        for i in range(number): # loop all number of points in the random samples
            print i
            x_index = round(samples_x[i])
            y_index = round(samples_y[i])
            z_index = round(samples_z[i])
            if image[x_index, y_index, z_index] > threshold:
                keep_points.append([samples_x[i], samples_y[i], samples_z[i]])
                
        return np.array(keep_points)
            
                  
    def GenSamples(self,rck_in,lim):
        """
        this methodology of sample generator randomly selects the points from surface points to generate poisson disk samples
        To control the sample growth direction, in hull function was applied to only keep the point within the convex hull
        
        """
        # the following lines describe the random selected point from the ISO surface points cloud
        # it orinally designed with the ISO surface point constrained generation
        # surf_index = random.randint(0,self.surf_points.shape[0]-1)
        # z = self.surf_points[surf_index][0]
        # x = self.surf_points[surf_index][1]
        # y = self.surf_points[surf_index][2]
        ## be careful about the arrangement of z, x, y, in turns.
        # pos = np.array([[z],[x],[y]])
        
        pos = np.array([[],[],[]])
        
        while(True):    
            z = random.randint(0,self.seg_arr.shape[0]-1)
            x = random.randint(0,self.seg_arr.shape[1]-1)
            y = random.randint(0,self.seg_arr.shape[2]-1)
            # the following line guarantes the first random point is in the middle of image histogram
            
            # the level is calculated firstly by the Otsu's method
            # in order to generate uniformly distributed points, and checker for active sample list must be added.
            if self.seg_arr[z][x][y] >= self.threshold:   
            #if self.seg_arr[z][x][y] >= self.threshold:
                pos = np.concatenate((pos,np.array([[z],[x],[y]])), axis = 1)
                if pos.shape[1] == 100:
                    break
                
        # print pos
        # initiate the sample list
        sample_array = np.array([[],[],[]])
        #sample_array = self.surf_points.T
        
        # initiate the active list
        active_array = np.array([[],[],[]])
        #active_array = self.surf_points.T
        
        # add the first random point into sample list
        #sample_array = np.append(sample_array,[[pos[0][0]],[pos[1][0]],[pos[2][0]]],axis=1)
        sample_array = np.concatenate((sample_array,pos),axis=1)
        # add the first random point into active list
        #active_array = np.append(active_array,[[pos[0][0]],[pos[1][0]],[pos[2][0]]],axis=1)
        active_array = np.concatenate((active_array,pos),axis=1)
        i = 0        
        Xlim = np.array([[0,self.seg_arr.shape[0]-1],[0,self.seg_arr.shape[1]-1],[0,self.seg_arr.shape[2]-1]])
                    # the sample_number controls how many samples generated at each query point
        if sample_key == 1:
            sample_number = 50
            
        if sample_key == 2:
            sample_number = 50
            lim = int(lim/20)
            print 'The sample number limits for python environment have been reduced to ' + str(lim)
                
        while len(active_array[0])!=0:
            tsample = np.transpose(sample_array)
            #tactive = np.transpose(active_array)       
            active_lotty = len(active_array[0])
            # calculating the length of active list
            query_number = np.random.choice(active_lotty)
            # generating a random number for the array
            query_point = np.transpose(np.array([[active_array[0][query_number]],[active_array[1][query_number]],[active_array[2][query_number]]]))
            xq = int(query_point[0][0])
            yq = int(query_point[0][1])
            zq = int(query_point[0][2])               
            
            # a distance function for controling the radius of PDS
            #rck = math.ceil(self.sizing[xq][yq][zq]) *  self.mesh_scale_factor
            # try no ceiling method
            #rck = self.sizing[xq][yq][zq] +  self.mesh_scale_factor
            
            rck = rck_in
            
            #if self.sizing[xq][yq][zq]>0 and self.sizing[xq][yq][zq]<rck_in:
            #    rck = int(rck_in/2)
            
            # a stepper function for controling the radius of PDS
            # self.seg_arr[xq][yq][zq]
            #if distance_value >0 and distance_value < 3  :
                 #rck = 3
            #elif segvalue >3 and segvalue <=6 :
                 #rck = 10
            #elif segvalue >6 :
                 #rck = 15   
            #else:
                 #rck = 10
            
            # looper control index and it controls the maximum number ran on the looper.                        
            i+=1
            # the following index controls the maximum 
            if i == lim:
                break
                
            ## Fortran Interface for distance check

            output = ds1.distance(tsample, query_point,rck,sample_number,3,Xlim)
            
            #print output[0]
            # this output is (N,3) dimension
            sz = output[1]  # the length sz == 50 as function input above
            # there are two ways to confine the sample generation, one is to check if there is no future samples generated,
            # another is to check the sample is within the distance function
            
            ###############  (1) check the if there are samples generated ######
            if sz == 0:
                new_active_array0 = np.delete(active_array[0],[query_number])
                new_active_array1 = np.delete(active_array[1],[query_number])
                new_active_array2 = np.delete(active_array[2],[query_number])
                del active_array
                active_array = np.array([new_active_array0,new_active_array1,new_active_array2])
            else:
            ###############  (2) check the samples generated ###################
                toutput = np.transpose(output[0]) # this is (3,N) dimension and it will be concatenated along the axis = 1
                # the following method checks the distance value for each sample
                toutput2 = self.dis_check(toutput,self.threshold)
                # print toutput2.shape[1]
                # the next is to check the generated sample if there is out of the distance function
                # if it is, then delete it from the active list.
                if toutput2.shape[1] != 0:
                    sample_array = np.append(sample_array,toutput2,axis=1)
                    active_array = np.append(active_array,toutput2,axis=1)
                else:
                    new_active_array0 = np.delete(active_array[0],[query_number])
                    new_active_array1 = np.delete(active_array[1],[query_number])
                    new_active_array2 = np.delete(active_array[2],[query_number])
                    del active_array
                    active_array = np.array([new_active_array0,new_active_array1,new_active_array2])
            
            self.active_size.append(active_array.shape[1])
            
            print "the progress: " + str(self.active_size[-1]) + " " + str(i)
            # try to control the sample growth rate
            #if int(self.active_size[-1]/10) - int(self.active_size[-2]) < 0:
                #break
            
            #sample_array = np.append(sample_array,[toutput[0][1:sz+1],toutput[1][1:sz+1],toutput[2][1:sz+1]],axis=1)
            #active_array = np.append(active_array,[toutput[0][1:sz+1],toutput[1][1:sz+1],toutput[2][1:sz+1]],axis=1)
            
            #xs, ys = toutput2.shape
            #check_list = self.in_hull(output[0], self.surf_points)
            #check_index = [item for item in range(1,len(check_list)) if check_list[item] == True] 
            ## the range is from 1 to 49, total 49 numbers, the first is itself, overlap item.
            #for j in xrange(len(check_index)):
                ## be careful about the dimension and shape of numpy array and usage of concatenate
                ## sample_array = np.array([[],[],[]]) this is the shape                
                #apoint_array = np.array([[toutput[0,check_index[j]]],[toutput[1,check_index[j]]],[toutput[2,check_index[j]]]])             
                #sample_array = np.append(sample_array,apoint_array,axis=1)
                #active_array = np.append(active_array,apoint_array,axis=1) 
                               
        self.samples = sample_array          
        
if __name__ == '__main__':  
    #print "This is a test."
    #from Visulization.DataVisu import DataVisulization               
    #import scipy.io as sio
    #import Visulization.DataVisu 
    ##reload(DataVisu)
    #from LoadDicom import ReadDicom     
    #from MeshSize import MeshSize
    
    #mat_contents = sio.loadmat('mesh_test.mat')
    #seg_arr1 = mat_contents['V']
    # load the image from dicom files
    # the following codes need to be optimized, since it uses the read dicom twice.
    
    #PathDicom1 = "/Users/junchaowei/Desktop/Python_DVC/CT_Data/AluminumBar_Ref_Def/Dicom_image1"
    #dicom_volume1 = ReadDicom(PathDicom1)
    #dicom_volume1.loadfiles()
    #image1 = dicom_volume1.DicArray 
    
    # added on
    #image_trans = MeshSize(image1)
    #image_trans.grey2binary(4000)
    #image_trans.distance_transform()   
    
    #Samples = DiskSamples(image1,image_trans.newimage,4000)
    
    #Samples.init_surface(10000)
    #Samples.GenSamples(20,20000)
    #point_arr = Samples.samples
    #DataVisu.DataVisulization(image1,4000).scatterplot(point_arr)
    
    #Samples.image_histogram(image1,16)
    #Samples.otsu_method()
    #print Samples.level
    import scipy.io as sio
    from scipy import interpolate
    import Visulization
    from scipy import ndimage
    from scipy.interpolate import griddata

    from mayavi import mlab
    import mayavi

    import SpatialSampling
    import Visulization
    import InputOutput
    import CorrelationFunction
    import Registration
    from scipy import spatial
    from scipy import linalg

    import time

    mat_contents1 = sio.loadmat('/Users/junchaowei/Desktop/Python_DVC2/UP_Research/WholeRegion/793-whole-LC-smoothed.mat')
    mat_contents2 = sio.loadmat('/Users/junchaowei/Desktop/Python_DVC2/UP_Research/WholeRegion/808-whole-LC-smoothed.mat')

    image1 = mat_contents1['volume'][250:1180,200:1180,185:260]
    image2 = mat_contents2['volume'][250:1180,200:1180,185:260]

    threshold1 = image1.mean()
    threshold2 = image2.mean()

    #threshold1 = 70
    #threshold2 = 70

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
        Samples2.GenSamples(40, 120000)
        point_arr2 = Samples2.samples
        sample_point2 = np.transpose(point_arr2)
        Visulization.DataVisulization(image1, threshold).scatterplot(point_arr2)
        return sample_point2

    dist_img1 = mesh_dist(image1, threshold1)
    point_arr1 = sampling_def(image1, dist_img1, threshold1)

