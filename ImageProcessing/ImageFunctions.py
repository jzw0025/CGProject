from PIL import Image
import numpy as np
from scipy import ndimage
from matplotlib import pyplot as plt
import os
import scipy.io as sio
from skimage.exposure import rescale_intensity
from skimage.filters import threshold_otsu, threshold_adaptive
from skimage.measure import EllipseModel
from skimage.segmentation import active_contour
import math
import matplotlib.pyplot as plt
from skimage import exposure

##############
# 
# must use the command in Mac OS to convert the tiff image stack into tiff images:
#  " convert multipage.tif single%03d.tif  "
#
##############

def readOCTstack(PathTiff): 
    
    if type(PathTiff) is not str:
        raise TypeError("the input address must be a string!")
        
    if not os.path.isdir(PathTiff):
        raise TypeError("the input string must be a valid directory!")

    lstFilesDCM = []  # create an empty list
    print "found tiff files:" + str(len(lstFilesDCM))
    for dirName, subdirList, fileList in os.walk(PathTiff,topdown=False):
        #print dirName, subdirList, fileList
        for filename in fileList:
            if ".tif" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName,filename))
                # load the first file for the sample
                print "read the file:" + str(filename)
                if len(lstFilesDCM) == 1:
                    #print os.path.join(dirName,filename)
                    #im = Image.open(os.path.join(dirName,filename))
                    single_im = plt.imread(os.path.join(dirName,filename))
                    print "read"
                    
    imx, imy = single_im.shape 
    print "the size of input image is:" + str(imx)+ "/" + str(imy)         
    ConstPixelDims = (imx, imy, len(lstFilesDCM))
    ArrayTiff = np.zeros(ConstPixelDims, dtype=float) # int array
    i=0          
    for filenameDCM in lstFilesDCM:  
        print "processing the slice #: " + str(i)
        i += 1
        #im_data = plt.imread(filenameDCM) 
        im = Image.open(filenameDCM)
        imarray = np.array(im)
        ArrayTiff[:, :, lstFilesDCM.index(filenameDCM)] = imarray # save the processed slice array into image

    return ArrayTiff
    
def readPartialOCTstack(PathTiff, start, end): 
    
    if type(PathTiff) is not str:
        raise TypeError("the input address must be a string!")
        
    if not os.path.isdir(PathTiff):
        raise TypeError("the input string must be a valid directory!")

    lstFilesDCM = []  # create an empty list
    print "found tiff files:" + str(len(lstFilesDCM))
    for dirName, subdirList, fileList in os.walk(PathTiff,topdown=False):
        #print dirName, subdirList, fileList
        for filename in fileList:
            if ".tif" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName,filename))
                # load the first file for the sample
                print "read the file:" + str(filename)
                if len(lstFilesDCM) == 1:
                    #print os.path.join(dirName,filename)
                    #im = Image.open(os.path.join(dirName,filename))
                    single_im = plt.imread(os.path.join(dirName,filename))
                    print "read"
                    
    imx, imy = single_im.shape 
    print "the size of input image is:" + str(imx)+ "/" + str(imy)   
          
    ConstPixelDims = (imx, imy, len(lstFilesDCM))
    ArrayTiff = np.zeros(ConstPixelDims, dtype=float) # int array
    Slice = lstFilesDCM[start:end+1] # reduce the slice vertically
    
    i=0          
    for filenameDCM in Slice:  
        print "processing the slice #: " + str(i)
        i += 1
        #im_data = plt.imread(filenameDCM) 
        im = Image.open(filenameDCM)
        imarray = np.array(im)
        ArrayTiff[:, :, lstFilesDCM.index(filenameDCM)] = imarray # save the processed slice array into image

    return ArrayTiff

def medOCT(array, size, local=None):
    """
    perform medium filter for the input array:
        array --- ndarray (int)
        size --- local medium filter kernel size (int)
        local --- address string 
    """
    plt.ioff()
    out_volume = np.empty_like(array)
    for i in range(array.shape[2]):
        print "processing the image slice: " + str(i)
        f = ndimage.filters.median_filter(array[:,:,i], size)
        if local and os.path.isdir(local): # save the images into the directory
            fig = plt.figure()
            plt.imshow(f, interpolation='nearest', cmap='Greys_r')
            fig.savefig(local+'/d'+ str(i)+'.tif')
            plt.close(fig)
        out_volume[:,:,i] = f
        
    return out_volume
    
def toBinary(volume, threshold):
    """
    volume --- 3D volume
    threshold --- float value
    """
    if type(volume) != np.ndarray:
        
        raise TypeError('the input must be numpy array!')
        
    binary1 = np.zeros_like(volume) # make a copy of array
    
    binary1[ volume > threshold ] = 1
    
    return binary1
    
def sliceThreshold(volume, block_size = 5):
    """
    convert slice into binary using adaptive local ostu method
    
    volume --- 3D volume
    block_size --- int value
    
    """
    if type(volume) != np.ndarray:
        
        raise TypeError('the input must be numpy array!')

    x, y, z = volume.shape
    
    segImg = np.empty_like(volume)
    
    for i in range(z):
        
        binary_adaptive = threshold_adaptive(volume[:,:,i], block_size, offset=0)
        
        segImg[:,:,i] = binary_adaptive
    
    return segImg

def DynamicRangeImage(Image, inRange, outRange=None, pad=False):
    """
    input:
        Image --- (N,d) numpy ndarray
        inRange --- (min, max) 
        outRange --- (min, max)
    """
    
    sx,sy,sz = Image.shape
    if pad:
        
        pad = 10

        output = np.empty((sx+2*pad, sy+2*pad, sz+2*pad)) 
        #print "the padded volume size is: " + str(output.shape)
        
        volume = np.empty_like(Image)
        
        for i in range(volume.shape[2]):
            
            volume[:,:,i] = rescale_intensity(Image[:,:,i], in_range=inRange)
    
        output[pad:-pad, pad:-pad, pad:-pad] = volume
        print "the input volume size is: " + str(Image.shape)
        print "the padded volume size is: " + str(output.shape)
        return output
        
    else:
        
        
        volume = np.empty_like(Image)
        
        for i in range(volume.shape[2]):
            
            volume[:,:,i] = rescale_intensity(Image[:,:,i], in_range=inRange)
    
        print "the input volume size is: " + str(Image.shape)
        print "the padded volume size is: " + str(volume.shape)
        return volume
        
def HistEqual(Image):
    
    volume = np.empty_like(Image)
    
    for i in range(Image.shape[2]):
        
        volume[:,:,i] = exposure.equalize_hist(Image[:,:,i])
        
    return volume
    
def AdaptiveHistEqual(Image, Contrast=0.005):
    
    sx, sy, sz = Image.shape
    
    volume = np.empty((sx, sy, sz),dtype=float)
    
    for i in range(Image.shape[2]):
        print i
        volume[:,:,i] = exposure.equalize_adapthist(Image[:,:,i],clip_limit=Contrast)
        
    return volume
    
def SpaceMedianFilter(Image, size):
    
    if size >7:
        print "the size needed to be less than 7, it will take a long time to run!"
    
    volume = np.empty_like(Image)
    
    volume = ndimage.filters.median_filter(Image, size)
        
    return volume
    
def ConvertImageSlice(volume, address=None):
    """
    volume --- ndarray 
    address --- folder address
    
    """
    if not address:
        print "no address is imported: please make a fold on desktop named: image"
        address = '/Users/junchaowei/Desktop/image/'
    
    for i in range(volume.shape[2]):
        print "converting slice:" + str(i)
        plt.imshow(volume[:,:,i], interpolation='nearest', cmap='Greys_r')
        plt.savefig(address+ str(i) +'.tif')
        
    print "finished converting image slices"
    
class LineBuilder:
    """
    interactive function for the line model
    """
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        print('click', event)
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0
    
def KeyPoints(array):
    """
    this return the Liner object with selected points in the 'array' image
    """
    fig = plt.figure()
    plt.imshow(array,cmap='Greys_r')  
    ax = fig.add_subplot(111)
    ax.set_title('click to build line segments')
    line, = ax.plot([0], [0])  # empty line
    linebuilder = LineBuilder(line)  
    plt.show()
    #x_cor = linebuilder.xs
    #y_cor = linebuilder.ys
    #print x_cor,y_cor
    return linebuilder
    
def FindEllipse(points, number=100):
    """
    this method applies the ellipse model to the input points
    points --- (N,2) inputs array
    number --- number of points in the estimated ellipse 
    -----
    out ---  return points on the ellipse
    """
    if type(points) is not np.ndarray:
        raise TypeError("the input must be a nd-array")
        
    if points.shape[1] !=2:
        print "the input array must be (N,2)" 
        return
    
    model = EllipseModel() # create an object
    model.estimate(points)
    out = model.predict_xy(np.linspace(0,2*math.pi,number))
    return out
    
def LogicRegion(array, points):
    """
    this function creates logic region using the geometric shapes;
    the points inside of the shape are mapped into "unit ones";
    the points outside of the shape remains zeros
    array ---- nd array
    points --- (N,2) geometric input, ellipse point or circle points
    ---
    out: Index_logic --- logical map with size equals to the array input
    """
    imx, imy = array.shape
    mx, my = np.meshgrid(np.linspace(0,imx-1,imx), np.linspace(0,imy-1,imy))
    # meshgrid change the orientation of x,y: x becomes the horizontal axis, y becomes the vertical
    # this change would affect the return statement of reshape.
    Cord_xy = np.vstack([mx.flatten(), my.flatten()]).T # column major vector
    Inhull_point = in_hull(Cord_xy, points)
    Index_logic = Inhull_point.reshape(imx, imy)
    #plt.imshow(Index_logic)
    #plt.show()
    return Index_logic    
  
def ActiveEllipse(array, points):
    """
    read the docs for the active_contour with tuned parameters
    http://scikit-image.org/docs/dev/auto_examples/edges/plot_active_contours.html
    """
    snake = active_contour(array, points, w_edge=-1, w_line=1, gamma=0.001, beta=10) 
    #plt.imshow(a)
    #plt.plot(snake[:,0],snake[:,1])
    return snake
    
def ImageRescale255(array):
    """
    this function recalculates the range of image(0.0,1.0) into (0,255)
    """
    maxa = array.max()
    mina = array.min()
    volume = 255.0*(array-mina)/(maxa-mina)
    return volume.astype(int)
    
def ImageRescale1(array):
    """
    this function recalculates the range of image (0,255) into (0.0,1.0)
    """
    maxa = array.max()
    mina = array.min()
    volume = 1.0*(array-mina)/(maxa-mina)
    return volume.astype(float)
    
def AspectRatioPoins(Points,scale):
    """
    rescale point (N,3) for 3D printing
    
    scale--- tuple (min, max)
    
    """
    
    min0 = Points[:,0].min()
    max0 = Points[:,0].max()
    
    min1 = Points[:,0].min()
    max1 = Points[:,0].max()
    
    min2 = Points[:,0].min()
    max2 = Points[:,0].max()
    
    normalized = np.empty_like(Points)
    normalized[:,0] = (Points[:,0]-min0)/(max0-min0)
    normalized[:,1] = (Points[:,1]-min1)/(max1-min1)
    normalized[:,2] = (Points[:,2]-min1)/(max1-min1)
    
    # resacale
    normalized[:,0] = normalized[:,0]*(scale[1]-scale[0])+scale[0]
    normalized[:,1] = normalized[:,1]*(scale[1]-scale[0])+scale[0]
    normalized[:,2] = normalized[:,2]*(scale[1]-scale[0])+scale[0]
    
    return normalized
    
def hist_match(source, template):
    """
    Adjust the pixel values of a grayscale image such that its histogram
    matches that of a target image

    Arguments:
    -----------
        source: np.ndarray
            Image to transform; the histogram is computed over the flattened
            array
        template: np.ndarray
            Template image; can have different dimensions to source
    Returns:
    -----------
        matched: np.ndarray
            The transformed output image
    """

    oldshape = source.shape
    source = source.ravel()
    template = template.ravel()

    # get the set of unique pixel values and their corresponding indices and
    # counts
    s_values, bin_idx, s_counts = np.unique(source, return_inverse=True,
                                            return_counts=True)
    t_values, t_counts = np.unique(template, return_counts=True)

    # take the cumsum of the counts and normalize by the number of pixels to
    # get the empirical cumulative distribution functions for the source and
    # template images (maps pixel value --> quantile)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]

    # interpolate linearly to find the pixel values in the template image
    # that correspond most closely to the quantiles in the source image
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)

    return interp_t_values[bin_idx].reshape(oldshape)
    
def IntensityNormalization(window):  
    #sx, sy, sz = window.shape
    #mean = window.mean()
    #stad = window.std()
    #outwindow = (window - mean)/stad
    maxIntensity = window.max()
    minIntensity = window.min()
    window = (window - minIntensity) / (maxIntensity - minIntensity)
    return window
    
def applyLogic(Image, logic):
    
    sx, sy, sz = Image.shape
    
    volume = np.empty_like(Image)
    for i in range(sz):
        volume[:,:,i] = Image[:,:,i] * logic
        
    return volume
    
    
if __name__ == "__main__":
    print "this is test file!"
    out = readOCTstack("/Users/junchaowei/Desktop/Python_DVC2/UP_Research/WholeRegionRealData/For Junchao/M6_OD_125_C-scan_full")                    
    out2 = readPartialOCTstack("/Users/junchaowei/Desktop/Python_DVC2/UP_Research/WholeRegionRealData/For Junchao/M6_OD_baseline_C-scan_full",350,400)                    
    medOCT(out,15, local="/Users/junchaowei/Desktop/test")
    data = np.array([linebuilder.xs[1:], linebuilder.ys[1:]]).T    
    points = FindEllipse(data)    
    a = plt.imread("/Users/junchaowei/Desktop/Python_DVC2/UP_Research/WholeRegionRealData/For Junchao/M6_OD_125_C-scan_full/400.tiff")
    ra = rescale_intensity(a, in_range=(60,255))
    
    mra = ndimage.filters.median_filter(ra, 3)
    
    LogicRegion(a,points) 
    plt.imshow(ra,cmap='Greys_r')
    plt.plot(out[:,0],out[:,1])

    
    
    

      