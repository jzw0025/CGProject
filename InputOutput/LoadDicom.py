# -*- coding: utf-8 -*-
import numpy
import os   # import the operation system function
os.getcwd() # get the current working directory
import pydicom # import the pydicom package

class ReadDicom:
    """ This class read the dicom files inside the folder
     load the dicom files"""
    def __init__(self,directory): 
        " initialize the variables "
        if type(directory)!=str:
            print "unknown directory format!"     
        self.PathDicom = directory
        self.DicArray = numpy.array([]) 
        
    def loadfiles(self):
        lstFilesDCM = []  # create an empty list
        for dirName, subdirList, fileList in os.walk(self.PathDicom):
            for filename in fileList:
                if ".dcm" in filename.lower():  # check whether the file's DICOM
                    lstFilesDCM.append(os.path.join(dirName,filename))           
        # If you check the MyHead folder you’ll see that the .dcm files are named MR000000.dcm, MR000001.dcm, etc.
        # Therefore, the walk function will return them in order since they’re sorted lexicographically by the OS
        # Get ref file            
        RefDs = pydicom.read_file(lstFilesDCM[0])
        # adding the attribute in order to avoid error message
        RefDs.SamplesPerPixel=1
        # ArrayDicom[:, :, 0] = RefDs.PixelData
        # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
        ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))
        ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]),float(RefDs.PixelSpacing[1]))
        # we simply use numpy.arange, ConstPixelDims, and ConstPixelSpacing to calculate axes for this array
        x = numpy.arange(0.0, (ConstPixelDims[0]+1)*ConstPixelSpacing[0], ConstPixelSpacing[0])
        y = numpy.arange(0.0, (ConstPixelDims[1]+1)*ConstPixelSpacing[1], ConstPixelSpacing[1])
        z = numpy.arange(0.0, (ConstPixelDims[2]+1)*ConstPixelSpacing[2], ConstPixelSpacing[2])
        # load the grey scale value from dicom files
        # The array is sized based on 'ConstPixelDims'
        ArrayDicom = numpy.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)
        # loop through all the DICOM files
        for filenameDCM in lstFilesDCM:
        # read the file
            RefDs = pydicom.read_file(filenameDCM)
            RefDs.SamplesPerPixel=1
            # store the raw image data
            ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = RefDs.pixel_array 
            self.DicArray = ArrayDicom
        
if __name__ == "__main__":
    print "this is a test."
    PathDicom2 = "/Users/junchaowei/Desktop/Python_DVC2/CT_Data/dicom512def"
    Test = ReadDicom(PathDicom2)
    Test.loadfiles()
        
        