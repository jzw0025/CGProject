import numpy as np
import scipy.ndimage.filters as filters
import scipy.io as sio
from scipy import ndimage
from math import *
from scipy import signal
import Visualization
import sys
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from mayavi import mlab
from mayavi.mlab import *

load = sio.loadmat('/Users/junchaowei/Desktop/Pack12122016/clean_000.mat')
image1 = load['par1'][200:300,200:300,50:100]

load = sio.loadmat('/Users/junchaowei/Desktop/Pack12122016/clean487_reg1.mat')
image2 = load['par1'][200:300,200:300,50:100]

########### create a test volume ###########
#test_image1 = np.zeros((50,50,50))
#rnd_pts = 50*np.random.random([25,3])
#rnd_pts2 = rnd_pts.astype(int)
#for i in range(rnd_pts.shape[0]):
#    test_image1[rnd_pts2[i][0]][rnd_pts2[i][1]][rnd_pts2[i][2]] = 1   
#new_test_image1 = ndimage.binary_dilation(test_image1,iterations=3).astype(test_image1.dtype)
#smooth_new_test_image1 = ndimage.gaussian_filter(new_test_image1,2)
#v1 = Visualization.DataVisulization(smooth_new_test_image1,0.2) # initial reference image
#v1.contour3d()
#image1 = smooth_new_test_image1
#plt.figure(1)
#plt.imshow(image1[:,:,25])
#plt.show()

sx, sy, sz = image1.shape
    
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
        
def createPoints(image1):
    sx, sy, sz = image1.shape
    global lx,ly,lz
    ### matlab fspeical function python port:
    ### references: http://stackoverflow.com/questions/17190649/how-to-obtain-a-gaussian-filter-in-python
        
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
    
    def findLocalMaximum(val, radius):
        """
        """
        def mask(radius): # create ball mask
            """
            this is the help function, it is a distance mask,on which the ones are placed.
            """
            y,x,z = np.ogrid[-radius: radius+1, -radius: radius+1, -radius: radius+1]
            mask = x**2 + y**2 +z**2 <= radius**2
            return mask
            
        val_x = int(val.shape[0])
        val_y = int(val.shape[1])
        val_z = int(val.shape[2])
        radius = int(radius)
        max_local   = np.zeros((val_x,val_y,val_z))
        val_enlarge = np.zeros((val_x+2*radius,val_y+2*radius,val_z+2*radius))
        val_mask    = np.zeros((val_x+2*radius,val_y+2*radius,val_z+2*radius))
    
        val_enlarge[radius:val_x+radius,
                           radius:val_y+radius,
                           radius:val_z+radius] = val
                    
        val_mask[radius:val_x+radius,
                        radius:val_y+radius,
                        radius:val_z+radius] = 1
                
        msk = mask(radius) 
    
        dimX = []
        dimY = []
        dimZ = []
        index=0
        for i in range(val_x):
            print "the processed kernal is:" + str(index*100.0/val_x/val_y/val_z)
            for j in range(val_y):
                for k in range(val_z):
                    index+=1
                    val_ref = val[i, j, k]
                    neigh_val = val_enlarge[i:i+2*radius+1, j:j+2*radius+1, k:k+2*radius+1] # X-A,X+A,Y-B,Y+B,Z-C,Z+C: i+radius-Radius:i+3*radius+1-Radius
                    neigh_mask = val_mask[i:i+2*radius+1, j:j+2*radius+1, k:k+2*radius+1]*msk
                    #neigh_mask = msk
                    #neigh_sort = neigh_val[neigh_mask==1].flatten()
                    #neigh_sort.sort()
                    #if val_ref==neigh_sort[-1] and val_ref>neigh_sort[-2]: # saving time
                    if val_ref == neigh_val[neigh_mask==1].max(): # if this value is the local maximum value.
                        dimX.append(i)
                        dimY.append(j)
                        dimZ.append(k)
                        max_local[i,j,k] = val_ref
                        
        return dimX,dimY,dimZ,max_local
                                  
    #################### harris ###################   
    # scale paramter
    sigma_begin = 1.5
    sigma_step = 1.2
    sigma_nb = 6
    sigma_array = []
    for i in range(sigma_nb):
        sigma_array.append(sigma_step**(i)*sigma_begin) 
        
    harris_pts = [[],[],[],[]]
    for i in range(len(sigma_array)):
        
        # scale (standard deviation)
        s_I = sigma_array[i] # integration scale
        s_D = 0.7*s_I # derivative scale 0.7
        
        # image derivations
        Ix = filters.gaussian_filter(image1, sigma=s_D, order=[1,0,0])
        Iy = filters.gaussian_filter(image1, sigma=s_D, order=[0,1,0])
        Iz = filters.gaussian_filter(image1, sigma=s_D, order=[0,0,1])
        
        size = max(1, floor(6*s_I+1))
        
        # auto-correlation matrix
        g = gauss3D(shape=(size,size,size),sigma=s_I) # set up 3D gaussian kernal 
        #Ix2 = filters.convolve(Ix*Ix,g)
        #Iy2 = filters.convolve(Iy*Iy,g)
        #Ixy = filters.convolve(Ix*Iy,g)
        Ix2 = signal.fftconvolve(Ix*Ix,g,mode='same')
        Iy2 = signal.fftconvolve(Iy*Iy,g,mode='same')
        Iz2 = signal.fftconvolve(Iz*Iz,g,mode='same')
        Ixy = signal.fftconvolve(Ix*Iy,g,mode='same')
        Ixz = signal.fftconvolve(Ix*Iz,g,mode='same')
        Iyz = signal.fftconvolve(Iy*Iz,g,mode='same')
        
        #Hxyz = np.array()
        # interest point response
        k = 0.06 
        #cim = (Ix2 + Iy2 +Iz2)**3/(Ix2*Iy2*Iz2 + Ixy*Iyz*Ixz + Ixy*Iyz*Ixz - Ixy*Ixy*Iz2 - Iyz*Iyz*Ix2 - Ixz*Ixz*Iy2) # 3 dimensional from 2D criterion is not working properly
        #cim = (Ix2*Iy2*Iz2 + Ixy*Iyz*Ixz + Ixy*Iyz*Ixz - Ixy*Ixy*Iz2 - Iyz*Iyz*Ix2 - Ixz*Ixz*Iy2) - k*(Ix2 + Iy2 +Iz2)**3 
        
        cim = Ix2*Iy2 - Ixy**2 - k*(Ix2 + Iy2)**2
        #cim2 = Ix2*Iz2 - Ixz**2 - k*(Ix2 + Iz2)**2
        #cim3 = Iy2*Iz2 - Iyz**2 - k*(Iy2 + Iz2)**2
        #cim = cim1 + cim2 +cim3
        
        ix,iy,iz,ival = findLocalMaximum(cim,3*s_I)
        
        # set the threshold 2% iif the maximum
        t = 0.2*ival.max()
        ix,iy,iz = np.where(ival>t)
        
        # build array of interest points
        n = ix.shape[0]
        
        harris_pts[0].extend(ix.tolist())
        harris_pts[1].extend(iy.tolist())
        harris_pts[2].extend(iz.tolist())
        harris_pts[3].extend([i]*n)
        print len(harris_pts[0])
                                    
    #################### Laplace ###################   
    # compute scale-normalized laplacian operator
    laplace_snlo = np.zeros((sx,sy,sz,len(sigma_array)))
    for i in range(len(sigma_array)):
        s_L = sigma_array[i] # scale
        LoG = lapGauss3D(shape=(floor(6*s_I+1),floor(6*s_I+1),floor(6*s_I+1)), std=s_L)
        laplace_snlo[:,:,:,i] = signal.fftconvolve(image1,LoG,mode='same')
        
    # verify for each of the initial points whether the LoG attains a maximum at the scale of the point
    nub = len(harris_pts[0])
    points = [[],[],[],[]] 
    for i in range(nub):
        ixx = harris_pts[0][i]
        iyy = harris_pts[1][i]
        izz = harris_pts[2][i]
        s = harris_pts[3][i]
        
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
                
    return points  
    
points1 = createPoints(image1) 
points2 = createPoints(image2) 
    
v1 = Visualization.DataVisulization(image1,image1.mean()) # initial reference image
v1.contour3d()

points3d(points1[0], points1[1], points1[2], points1[3], colormap="Greens", scale_factor=1)
points3d(points2[0], points2[1], points2[2], points2[3], colormap="pink", scale_factor=1)

def getDescriptor(points):
    ### icosahedron construction ###
    phi = (1.0 + sqrt(5.0))/2
    vertices =[[-1, phi,0],
                    [1, phi, 0],
                    [-1,-phi,0],
                    [1, -phi,0],
                    [0,-1, phi],
                    [0, 1, phi],
                    [0,-1,-phi],
                    [0, 1,-phi],
                    [phi,0, -1],
                    [phi, 0, 1],
                    [-phi,0,-1],
                    [-phi,0, 1]]
    
    elements = [[0,11,5],
                        [0,5, 1],
                        [0,1, 7],
                        [0,7,10],
                        [0,10,11],
                        [1,5, 9],
                        [5,11,4],
                        [11,10,2],
                        [10,7,6],
                        [7,1,8],
                        [3,9,4],
                        [3,4,2],
                        [3,2,6],
                        [3,6,8],
                        [3,8,9],
                        [4,9,5],
                        [2,4,11],
                        [6,2,10],
                        [8,6,7],
                        [9,8,1]]
                
    vertices = np.array(vertices)
    elements = np.array(elements)            
    center = []
    for i in range(len(elements)):   
        center.append((vertices[elements[i][0]] + vertices[elements[i][1]] + vertices[elements[i][2]])/3.0)
    center = np.array(center)
    
    normals = []
    for i in range(center.shape[0]):
        norm = sqrt(sum(center[i,:]**2)) 
        normals.append(center[i,:]/norm)   
    normals = np.array(normals) # 20 bins
    
    #fig = mlab.figure(2)
    #mlab.clf()
    #mesh = mlab.triangular_mesh(vertices[:,0], vertices[:,1], vertices[:,2], elements)
    #cursor3d = mlab.points3d(0., 0., 0., mode='axes',
    #                                     color=(0, 0, 0),
    #                                     scale_factor=0.5)
    #points3d(center[:,0], center[:,1], center[:,2], colormap="Greens", scale_factor=0.1)
    
    #### bins and histogram
    descriptor = np.empty((len(points[0]),20*8))
    for n_points in range(len(points[0])): 
        x = points[0][n_points]
        y = points[1][n_points]
        z = points[2][n_points]
        s = floor(points[3][n_points])
        
        size = int(floor(6*s+1))  
        radius = int(size-1)/2 # half width of size
        # auto-correlation matrix
        weights =  gauss3D(shape=(size,size,size),sigma=(size-1)/4.0) # set up 3D gaussian kernal 
        
        pIx = np.empty((sx+2*radius,sy+2*radius,sz+2*radius))
        pIx[radius:-radius,radius:-radius,radius:-radius] = Ix
        
        pIy = np.empty((sx+2*radius,sy+2*radius,sz+2*radius))
        pIy[radius:-radius,radius:-radius,radius:-radius] = Iy
        
        pIz = np.empty((sx+2*radius,sy+2*radius,sz+2*radius))
        pIz[radius:-radius,radius:-radius,radius:-radius] = Iz
        
        wIx = pIx[x:x+2*radius+1,y:y+2*radius+1,z:z+2*radius+1]*weights # -radius + radius = 0
        wIy = pIy[x:x+2*radius+1,y:y+2*radius+1,z:z+2*radius+1]*weights  # +radius + radius = 2*radius
        wIz = pIz[x:x+2*radius+1,y:y+2*radius+1,z:z+2*radius+1]*weights 
        
        quand = [[0,radius,0,radius,0,radius],
                      [0,radius,0,radius,radius,size],
                      [0,radius,radius,size,0,radius],
                      [0,radius,radius,size,radius,size],
                      [radius,size,0,radius,0,radius],
                      [radius,size,0,radius,radius,size],
                      [radius,size,radius,size,0,radius],
                      [radius,size,radius,size,radius,size]] # Eight Quandrant
        
        bins = np.zeros((8,20)) # this arrangement is the same as normals
        q_index = 0
        for q in range(len(quand)):
            for i in range(quand[q][0],quand[q][1]):
                for j in range(quand[q][2],quand[q][3]):
                    for k in range(quand[q][4],quand[q][5]):
                        if wIx[i,j,k] and wIy[i,j,k] and wIz[i,j,k]:
                            mag = sqrt(wIx[i,j,k]**2+wIy[i,j,k]**2+wIz[i,j,k]**2)
                            norm = sqrt(sum(np.array([wIx[i,j,k],wIy[i,j,k],wIz[i,j,k]])**2))
                            dir_vector = np.array([wIx[i,j,k],wIy[i,j,k],wIz[i,j,k]])/norm
                            cosd_mag = np.dot(dir_vector,normals.T)
                            #print sum(np.dot(dir_vector,normals.T)>0.74535599249992979) 
                            index = cosd_mag > 0.74535599249992979 # 0.74535599249992979 is the angle of nearest norms
                            bins[q_index,:] += mag*cosd_mag*index
            q_index+=1
        descriptor[n_points,:] = bins.flatten()
        
    return descriptor

#for inor in range(descriptor.shape[0]):
#    normd = sqrt(sum(descriptor[inor,:]**2))
#    descriptor[inor,:] = descriptor[inor,:]/normd
#                                       
#for i in range(8):
#    plt.figure(2)
#    plt.plot(range(160),descriptor[10,:])


d1 = getDescriptor(points1)
d2 = getDescriptor(points2)

index_target = []
index_source = []
for pi in range(d1.shape[0]):
    t_p1 = d1[pi,:]
    ds = +inf
    index_target.append(pi)
    for ti in range(d2.shape[0]):
        temp = sqrt(sum((t_p1-d2[ti,:])**2))
        if temp < ds:
            ds = temp
            target_index = ti
    index_source.append(target_index) # must be at outside of loop 
    
v1 = Visualization.DataVisulization(image2,image2.mean()) # initial reference image
v1.contour3d()
  
# points1[0] --- x axis
# points1[1] --- y axis
# points1[2] --- z axis


points3d(points1[0], points1[1], points1[2], points1[3], colormap="Greens", scale_factor=0.5)
points3d(points2[0], points2[1], points2[2], points2[3], colormap="pink", scale_factor=0.5)

for ri in range(len(index_target)):
    x, y, z = [], [], []
    x.append(points1[0][index_target[ri]])
    x.append(points2[0][index_source[ri]])

    y.append(points1[1][index_target[ri]])
    y.append(points2[1][index_source[ri]])
    
    z.append(points1[2][index_target[ri]])
    z.append(points2[2][index_source[ri]])
    
    x = np.array(x).astype(float)
    y = np.array(y).astype(float)
    z = np.array(z).astype(float)
    plot3d(x, y, z, tube_radius=0.25, colormap='Spectral')


#x.extend([points1[0][0],points2[0][0]])
#y.extend([points1[1][0],points2[1][0]])
#z.extend([points1[2][0],points2[2][0]])

