import numpy
import os   # import the operation system function
os.getcwd() # get the current working directory
#from PIL import Image
import matplotlib.pyplot as plt

class ReadTiff:
    """ This class read the dicom files inside the folder
     load the dicom files"""
    def __init__(self,directory): 
        " initialize the variables "
        if type(directory)!=str:
            print "unknown directory format!"     
        self.PathDicom = directory
        self.single_im = []
        self.ArrayTiff = []

    def loadfiles(self):
        lstFilesDCM = []  # create an empty list
        for dirName, subdirList, fileList in os.walk(self.PathDicom):
            #print dirName, subdirList, fileList
            for filename in fileList:
                if ".tif" in filename.lower():  # check whether the file's DICOM
                    lstFilesDCM.append(os.path.join(dirName,filename))
                    # load the first file for the sample
                    if len(lstFilesDCM) == 1:
                        #print os.path.join(dirName,filename)
                        #im = Image.open(os.path.join(dirName,filename))
                        self.single_im = plt.imread(os.path.join(dirName,filename))
                        #print self.single_im
        
        imx, imy = self.single_im.shape 
        print imx,imy          
        ConstPixelDims = (imx, imy, len(lstFilesDCM))
        self.ArrayTiff = numpy.zeros(ConstPixelDims, dtype=self.single_im.dtype)
        for filenameDCM in lstFilesDCM:
        # read the file
            #im_data = Image.open(filenameDCM)
            im_data = plt.imread(filenameDCM)
            # store the raw image data
            self.ArrayTiff[:, :, lstFilesDCM.index(filenameDCM)] = numpy.array(im_data)
            
if __name__ == '__main__':
    from DataVisu import DataVisulization
    from PoissonDisk import DiskSamples
    from scipy import ndimage
    #import LoadTiff as LT
    from DelaunayMesh import MeshDelaunay
    from DelaunayMesh import MeshOptimizer
    
    # the following picture is 16 bit
    #tiff_dir = "/Users/jerry/Desktop/1533_5mmHg_700tdTomato_1-24um_zstack"
    
    # new image
    tiff_dir = "/Users/junchaowei/Desktop/Python_DVC/CT_Data/1533_5mmHg_800SHG_1-24um_zstack"
    # the following picture is also 16 bit 
    #tiff_dir = "/Users/jerry/Desktop/1533_5mmHg_800SHG_1-24um_zstack"
    Tiffs = ReadTiff(tiff_dir)
    Tiffs.loadfiles()
    print Tiffs.ArrayTiff
    #Samples = DiskSamples(Tiffs.ArrayTiff,[],[])
    #Samples.image_histogram(Tiffs.ArrayTiff,16)
    #Samples.otsu_method()
    #print Samples.level
    #Tiffs.ArrayTiff
    med_denoised = ndimage.gaussian_filter(Tiffs.ArrayTiff,3)
    # el = ndimage.generate_binary_structure(3, 1) # structing element 3D
    newimage = numpy.zeros_like(med_denoised)
    # the follwing converts the zeros array newimage into a binary array
    
    threshold = 6000 # choosing the separated level
    
    newimage[med_denoised < threshold ] = 0
    newimage[med_denoised >= threshold] = 1
    
    # perform multiple image openings
    ndimage.binary_opening(newimage, structure=numpy.ones((3,3,3))).astype(numpy.int)
    #ndimage.binary_opening(newimage, structure=numpy.ones((3,3,3))).astype(numpy.int)
    #ndimage.binary_opening(newimage, structure=numpy.ones((3,3,3))).astype(numpy.int)
    
    plot = DataVisulization(newimage,0.5) # 0.5 is choosen as to visulize binary image
    plot.contour3d()
    
    # the following line consumes a lot of computation time@
    newimage2 = ndimage.distance_transform_edt(newimage)
    
    Samples = DiskSamples(Tiffs.ArrayTiff,newimage2,6000)
    Samples.GenSamples(8,40000)
    point_arr = Samples.samples
    DataVisulization(newimage,25552).scatterplot(point_arr)
    
    mesh2 = MeshDelaunay(point_arr)
    mesh2.mesh()
    mesh2.alpha_shape(8)
    mesh2.mesh_triangle()
    
    Mesh_optimizer = MeshOptimizer(point_arr,mesh2.ntri)
    Mesh_optimizer.edge_connector_smooth()
    
    mesh2.view(mesh2.polydata(Mesh_optimizer.points.T))  
    
    #import matplotlib
    #import matplotlib.pyplot as plt

    #from skimage import data
    #from skimage.morphology import disk
    #from skimage.filters import threshold_otsu, rank
    #from skimage.util import img_as_ubyte
    
    #radius = 300
    #selem = disk(radius)
    #local_otsu = rank.otsu(Tiffs.ArrayTiff[:,:,50], selem)
    
                    
                    