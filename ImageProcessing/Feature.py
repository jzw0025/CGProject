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

def getDescriptor(image1, points):
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
    s = floor(0.7*points[3][0]/3) # this has been controlled by the harris points generator
    sx, sy, sz = image1.shape
    
    Ix = filters.gaussian_filter(image1, sigma=s, order=[1,0,0]) # save this preiously can make a faster computation
    Iy = filters.gaussian_filter(image1, sigma=s, order=[0,1,0])
    Iz = filters.gaussian_filter(image1, sigma=s, order=[0,0,1])
        
    descriptor = np.empty((len(points[0]),20*8))
    for n_points in range(len(points[0])): 
        print n_points
        x = points[0][n_points]
        y = points[1][n_points]
        z = points[2][n_points]
        
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
