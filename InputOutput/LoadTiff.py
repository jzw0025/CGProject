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
    from DelaunayMesh import MeshDelaunay
    from DelaunayMesh import MeshOptimizer
    
    # the following picture is 16 bit    
    # new image
    tiff_dir = " directory"
    # the following picture is also 16 bit 
    Tiffs = ReadTiff(tiff_dir)
    Tiffs.loadfiles()
    print Tiffs.ArrayTiff
                    