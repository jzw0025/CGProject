import numpy as np
from scipy import ndimage
# itertools is for using the looper three indexes i,j,k
# import itertools

"""
   This module mainly focus on creating a mesh size that adapted to the greyscale value,
   generatly, the seeds are according to the distance from ISO surface.
"""
class MeshSize:
    """
    this class creates the mesh seeds density for the usage of PoissonDisk Module
    
    """
    def __init__(self,ref_image):
        self.image = ref_image
        self.newimage = []
        
    def grey2binary(self,threshold):
        """
        this method converts the greyscale image into a binary image using the threshold method
        
        """
        ######## fast method #########
        # optinal 1:
        # image2 = np.where(self.image < threshold, self.image, 1)  
        # this function keep the elements less than the threshold, and replace the 
        # element greater than the threshod with 1
        
        # note: if it uses the self.image[self.image < threshold] instead of the self.newimage[]...
        # the original data array self.image will be modified, thus it will need double importing procedure for producing self.image
        
        self.newimage = np.zeros_like(self.image)
        
        # the follwing converts the zeros array newimage into a binary array
        self.newimage[self.image < threshold ] = 0
        self.newimage[self.image >= threshold] = 1
        
        ######## origin slow method(looper) #########
        #x_size, y_size, z_size = self.image.shape
        #self.newimage = np.zeros_like(self.image)
        #for x, y, z in itertools.product(xrange(x_size),xrange(y_size), xrange(z_size)):
            #if self.image[x][y][z] >= threshold:
                #self.newimage[x][y][z] = 1
                
    def region2binary(self, threshold1, threshold2):
        """
        this method converts only regional greyscale value into a binary format 
        """
        self.newimage = np.zeros_like(self.image)
        self.newimage[self.image < threshold1] = 0
        self.newimage[self.image > threshold2] = 0
        self.newimage[np.logical_and(self.image > threshold1,self.image < threshold2)] = 1
                
    def distance_transform(self):
        """
        this method converts the binary image into a distance image
        
        """
        image2 = self.newimage
        # the distance_transform_edit function transform the binary image into a distrance form image
        # for the detail: http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.morphology.distance_transform_edt.html
        self.newimage = ndimage.distance_transform_edt(image2)
        
if __name__ == '__main__':   
    from LoadDicom import ReadDicom
    PathDicom1 = "/Users/jerry/Desktop/Python_DVC/Dicom_image1"
    dicom_volume1 = ReadDicom(PathDicom1)
    dicom_volume1.loadfiles()
    image1 = dicom_volume1.DicArray 
    # the following two expressions would modify the image1 array into binary array
    image_trans = MeshSize(image1)
    #image_trans.grey2binary(4000)
    image_trans.region2binary(4000,6000)
    
    image_trans.distance_transform()
    print image_trans.image.max()
    print image_trans.newimage.max()
    print image_trans.newimage.min()
    print image_trans.newimage.mean()

        
        
                
            
            
        
            
            
            
