from PIL import Image
import numpy as np
from ImageProcessing import mediumGaussian

def readStack(dir, filted = None):
    """
    dir --- directory of image
    
    This method takes input directory of image, and read the image as a numpy stack
    """
    
    img = Image.open(dir)

    maxiter = 1000 # large default number for searching all the frames in tiff stack
    
    sx, sy = np.array(img).shape # all the later frame should have the same size as the first frame

    ImageData = []
    
    for i in range(maxiter):
        print "read image slice: " + str(i)
        try:
            img.seek(i)
             
            if filted is not None:
                imageSlice = mediumGaussian(np.array(img))
            else:
                imageSlice = np.array(img)
            
            ImageData.append(imageSlice) # frist store the data in the list, then stacking is faster than stacking array each time.
        
        except EOFError:
            # Not enough frames in img
            break
    
    return np.stack(ImageData, axis=2)
    

