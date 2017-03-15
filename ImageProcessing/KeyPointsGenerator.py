import numpy as np
import sys
import scipy.ndimage.filters as filters
from scipy import signal
from math import *
from scipy.ndimage import filters
from scipy import ndimage
from mayavi.mlab import *

def gauss3D(shape=(3,3,3),sigma=0.5):
        """
        3D gaussian mask - should give the same result as MATLAB's
        fspecial('gaussian',[shape],[sigma])
        the variances are equivalent in each direction
        """
        m,n,k = [(ss-1.)/2. for ss in shape]
        y,x,k = np.ogrid[-m:m+1,-n:n+1, -k:k+1]
        h = np.exp( -(x*x + y*y +k*k) / (2.*sigma*sigma) )
        h[ h < np.finfo(h.dtype).eps*h.max() ] = 0
        sumh = h.sum()
        if sumh != 0:
            h /= sumh
        return h

def lapGauss3D(shape=(3,3,3),std=0.5):
    """
    this function calculates the Laplacian of Gaussian
    h is the exponential kernal, h1 is the full expression of LoG
    """
    m,n,k = [(ss-1.)/2. for ss in shape]
    x = np.linspace(-m, m, 2*m+1)
    y = np.linspace(-n, n, 2*n+1)
    z = np.linspace(-k, k, 2*k+1)
    x, y, z = np.meshgrid(x, y, z)
    arg = -(x**2 + y**2 + z**2) / (2*std**2)
    h = np.exp(arg)
    h[h < sys.float_info.epsilon * h.max()] = 0
    h = h/h.sum() if h.sum() != 0 else h
    h1 = h*(x**2 + y**2 + z**2 - 2*std**2) / (std**4)
    return h1 - h1.mean()
    
def mask(radius): # create ball mask
    """
    this is the help function, it is a distance mask,on which the ones are placed.
    """
    y,x,z = np.ogrid[-radius: radius+1, -radius: radius+1, -radius: radius+1]
    mask = x**2 + y**2 +z**2 <= radius**2
    return mask

def findLocalMaximum(val, radius):
        """
        return the local maximum point
        val --- featured multidimension image
        radius --- local radius
        """
        radius = int(radius)
        
        val_x = int(val.shape[0])
        val_y = int(val.shape[1])
        val_z = int(val.shape[2])
    
        max_local   = np.zeros((val_x,val_y,val_z))
        val_enlarge = np.zeros((val_x+2*radius,val_y+2*radius,val_z+2*radius))
        val_mask    = np.zeros((val_x+2*radius,val_y+2*radius,val_z+2*radius))
        
        val_enlarge[radius:val_x+radius, radius:val_y+radius, radius:val_z+radius] = val # making artificial boundary    
        val_mask[radius:val_x+radius, radius:val_y+radius, radius:val_z+radius] = 1 # making artificial mask boundary          
        msk = mask(radius) 

        dimX, dimY, dimZ = [], [], []
        index=0
        for i in range(val_x):
            print "the processed kernal is:" + str(index*100.0/val_x/val_y/val_z)
            for j in range(val_y):
                for k in range(val_z):
                    index+=1
                    val_ref = val[i, j, k]
                    neigh_val = val_enlarge[i:i+2*radius+1, j:j+2*radius+1, k:k+2*radius+1] 
                    neigh_mask = val_mask[i:i+2*radius+1, j:j+2*radius+1, k:k+2*radius+1]*msk
                    if val_ref == neigh_val[neigh_mask==1].max(): # if this value is the local maximum value.
                        dimX.append(i)
                        dimY.append(j)
                        dimZ.append(k)
                        max_local[i,j,k] = val_ref              
        return dimX,dimY,dimZ,max_local

def createPoints(image1, density=0.1 ):
    sx, sy, sz = image1.shape
    # scale paramter
    sigma_begin = 1.5
    sigma_step = 1.2
    sigma_nb = 1
    sigma_array = []
    for i in range(sigma_nb):
        sigma_array.append(sigma_step**(i)*sigma_begin) 
    harris_pts = [[],[],[],[]]
    for i in range(len(sigma_array)):
        s_I = sigma_array[i] # integration scale
        s_D = 0.7*s_I # derivative scale 0.7
        Ix = filters.gaussian_filter(image1, sigma=s_D, order=[1,0,0])
        Iy = filters.gaussian_filter(image1, sigma=s_D, order=[0,1,0])
        #Iz = filters.gaussian_filter(image1, sigma=s_D, order=[0,0,1])
        size = max(1, floor(6*s_I+1))
        g = gauss3D(shape=(size,size,size), sigma=s_I) # set up 3D gaussian kernal 
        Ix2 = signal.fftconvolve(Ix*Ix,g,mode='same')        
        Iy2 = signal.fftconvolve(Iy*Iy,g,mode='same')
        #Iz2 = signal.fftconvolve(Iz*Iz,g,mode='same')
        Ixy = signal.fftconvolve(Ix*Iy,g,mode='same')
        #Ixz = signal.fftconvolve(Ix*Iz,g,mode='same')
        #Iyz = signal.fftconvolve(Iy*Iz,g,mode='same')
        k = 0.06 
        #cim = s_D**3*(Ix2*Iy2*Iz2 + Ixy*Iyz*Ixz + Ixy*Iyz*Ixz - Ixy*Ixy*Iz2 - Iyz*Iyz*Ix2 - Ixz*Ixz*Iy2) - k*(Ix2 + Iy2 +Iz2)**3 
        cim = Ix2*Iy2 - Ixy**2 - k*(Ix2 + Iy2)**2
        #cim2 = Ix2*Iz2 - Ixz**2 - k*(Ix2 + Iz2)**2
        #cim3 = Iy2*Iz2 - Iyz**2 - k*(Iy2 + Iz2)**2
        ix,iy,iz,ival = findLocalMaximum(cim,3*s_I)
        # set the threshold 2% iif the maximum
        t = density*ival.max()
        ix,iy,iz = np.where(ival>t)
        n = ix.shape[0] 
        harris_pts[0].extend(ix.tolist())
        harris_pts[1].extend(iy.tolist())
        harris_pts[2].extend(iz.tolist())
        harris_pts[3].extend([i]*n)
        print  "generator points number: " + str(len(harris_pts[0]))
                                 
    laplace_snlo = np.zeros((sx,sy,sz,len(sigma_array)))
    for i in range(len(sigma_array)):
        s_L = sigma_array[i] 
        LoG = lapGauss3D(shape=(floor(6*s_I+1),floor(6*s_I+1),floor(6*s_I+1)), std=s_L)
        laplace_snlo[:,:,:,i] = signal.fftconvolve(image1,LoG,mode='same')
    
    nub = len(harris_pts[0])
    points = [[],[],[],[]] 
    for i in range(nub):
        ixx = harris_pts[0][i]
        iyy = harris_pts[1][i]
        izz = harris_pts[2][i]
        s = harris_pts[3][i]
        if len(sigma_array) >=2:
            val = laplace_snlo[ixx,iyy,izz,s]
            if s>0 and s< len(sigma_array)-1:
                if val>laplace_snlo[ixx,iyy,izz,s-1] and val>laplace_snlo[ixx,iyy,izz,s+1]:
                    print i
                    points[0].append(ixx)
                    points[1].append(iyy)
                    points[2].append(izz)
                    points[3].append(3*sigma_array[s])
            elif s==0:
                if val>laplace_snlo[ixx,iyy,izz,1]:
                    points[0].append(ixx)
                    points[1].append(iyy)
                    points[2].append(izz)
                    points[3].append(3*sigma_array[s])
                    
            elif s==len(sigma_array)-1:
                if val>laplace_snlo[ixx,iyy,izz,len(sigma_array)-2]:
                    points[0].append(ixx)
                    points[1].append(iyy)
                    points[2].append(izz)
                    points[3].append(3*sigma_array[s]) 
        else:
                    points[0].append(ixx)
                    points[1].append(iyy)
                    points[2].append(izz)
                    points[3].append(3*sigma_array[s]) 
                
    return points 
    
def createPoints_sober(image1):
        Ix = ndimage.sobel(image1,0) 
        Iy = ndimage.sobel(image1,1)
        Iz = ndimage.sobel(image1,2)
        Ix2 = Ix*Ix
        Iy2 = Iy*Iy
        Iz2 = Iz*Iz
        Ixy = Ix*Iy
        Ixz = Ix*Iz
        Iyz = Iy*Iz
        k = 0.06 
        I_feature = (Ix2*Iy2*Iz2 + Ixy*Iyz*Ixz + Ixy*Iyz*Ixz - Ixy*Ixy*Iz2 - Iyz*Iyz*Ix2 - Ixz*Ixz*Iy2) - k*(Ix2 + Iy2 +Iz2)**3 
        return I_feature
        
def imagePading(image,padsize = 1):
    """
    inputs:
        image --- Nd array
    outputs:
        imageP  --- Nd array
    """
    sx,sy,sz = image.shape
    imageP = np.empty((sx+2*padsize, sy+2*padsize, sz+2*padsize))
    imageP[1:-1, 1:-1, 1:-1] = image
    
    return imageP
        
if __name__ == "__main__":
    print "this is a unit test!"
    import scipy.io as sio
    import matplotlib.pyplot as plt
    import skimage13 as skimage
    from skimage.transform import pyramid_gaussian
        
    out = createPoints(image1[100:300,100:300,0:100])
    out2 = createPoints(image2[100:300,100:300,0:100])

    x, y, z = out[0], out[1], out[2]
    points3d(x, y, z, colormap="Greens",scale_factor=10)
    import Visualization
    v1 = Visualization.DataVisulization(imagePading(image1[100:300,100:300,0:100]),0.2) # initial reference image
    v1.contour3d()
    
    x2, y2, z2 = out2[0], out2[1], out2[2]
    points3d(x2, y2, z2, colormap="flag", scale_factor=10)
    import Visualization
    v1 = Visualization.DataVisulization(imagePading(image2[100:300,100:300,0:100]),0.2) # initial reference image
    v1.contour3d()
    
    d1 = getDescriptor(image1[100:300,100:300,0:100],out)
    d2 = getDescriptor(image2[100:300,100:300,0:100],out2)
    
