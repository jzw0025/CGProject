import numpy as np
#import DccFort as FDCC
from numpy import unravel_index
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

import scipy.io as sio

import scipy.interpolate as interp
from scipy.interpolate import RegularGridInterpolator
import Visualization
import scipy.ndimage as ndimage

from pylab import imshow, show, get_cmap
from numpy import random

from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import scipy.io as sio
from scipy import linalg
from numpy.linalg import inv
from scipy import interpolate


class Correlation:
    "This class construct correlations"
    def __init__(self,point,seg_arr1,seg_arr2,subsize,subsearch):
        # Data Structure
        self.point_xyz = point
        self.subsize_fixed = subsize
        self.subsize_search = subsearch
        self.image1 = self.boundarypadding(seg_arr1,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
        self.image2 = self.boundarypadding(seg_arr2,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
        
        self.zdim = self.image1.shape[0]
        self.xdim = self.image1.shape[1]
        self.ydim = self.image1.shape[2]
        
        self.correlationz = np.zeros(shape=(self.zdim,self.xdim,self.ydim))
        self.correlationx = np.zeros(shape=(self.zdim,self.xdim,self.ydim))
        self.correlationy = np.zeros(shape=(self.zdim,self.xdim,self.ydim))

        self.point_z_dis = np.zeros(shape = (self.point_xyz.shape[0]))
        self.point_x_dis = np.zeros(shape = (self.point_xyz.shape[0]))
        self.point_y_dis = np.zeros(shape = (self.point_xyz.shape[0]))
        
        self.extrapo_z = []
        self.extrapo_x = []
        self.extrapo_y = []
        
        self._oridic = {}
        self._defdic = {}
        
        self.progress = 0.0
        self._sr = 1
        self.nx = np.linspace(-1,1,3)
        self.ny = np.linspace(-1,1,3)
        self.nz = np.linspace(-1,1,3)
        
        sx, sy, sz = self.image1.shape 
        x = np.arange(0, self.image1.shape[0])
        y = np.arange(0, self.image1.shape[1])
        z = np.arange(0, self.image1.shape[2])
        
        self.fff = RegularGridInterpolator((x, y, z), self.image1) # interpolating the original image1
        
    def DCC(self,subset1,subset2):
        Z, X, Y = subset1.shape
        Zf, Xf, Yf = subset2.shape
        # get the dimension of two subsets
        Z2 = Z // 2
        X2 = X // 2
        Y2 = Y // 2
        npad = ((Z2,Z2),(X2,X2),(Y2,Y2))
        psubset2 = np.pad(subset2,pad_width=npad, mode='constant', constant_values=0)
        Zfp, Xfp, Yfp = psubset2.shape
        ## size after zero-padding
        mean_subset1 = np.mean(subset1)
        #print subset1
        if mean_subset1 == 0.0:
            return (0,0,0)
        # mean of the subset1 (mean of flatten array)
        ####################  Fortran Interface  ##################
        forrij = FDCC.dcc(subset1,psubset2,Z,Zfp,Zf) 
        ###################  Python Interface  ###################  
        #V1 = subset1 - mean_subset1
        #V2 = V1**2 
        #rij = np.zeros(subset2.shape)
        #for z in range(Z2,Zfp-Z2):
            #for x in range(X2,Xfp-X2):
                #for y in range(Y2,Yfp-Y2):
                    ##csubset2z = psubset2.take(range(subset_zl,subset_zu + 1), axis=0) 
                    ##csubset2x = csubset2z.take(range(subset_xl,subset_xu + 1), axis=1)
                    ##csubset2y = csubset2x.take(range(subset_yl,subset_yu + 1), axis=2)
                    ##print subset_zl,subset_zu + 1,subset_xl,subset_xu + 1
                    #csubset2y = psubset2[z-Z2:z+Z2+1, x-X2:x+X2+1, y-Y2:y+Y2+1]
                    #mean_csubset2y = np.mean(csubset2y)                
                    #G1 = csubset2y - mean_csubset2y
                    #G2 = G1**2
                    #V2S = np.sum(V2)
                    #G2S = np.sum(G2)
                    #if V2S == 0 or G2S == 0:
                        ## just in case not dividing zeros
                        #corrivalue = 2.0
                    #else:                   
                        #corrivalue = 1.0 - np.sum(V1*G1)/np.sqrt(V2S)/np.sqrt(G2S)
                        #print corrivalue
                    #rij[z-Z2][x-X2][y-Y2] = corrivalue
                    ##print rij[z-Z2][x-X2][y-Y2]
        #print rij.min()
        #print forrij-rij

        R_index = unravel_index(forrij.argmin(),forrij.shape)          
        #Rinter = forrij[R_index[0]-self._sr : R_index[0]+self._sr+1,R_index[1]-self._sr : R_index[1]+self._sr+1,R_index[2]-self._sr : R_index[2]+self._sr+1]  
        corre_center = ([(Zf-1)/2,(Xf-1)/2,(Yf-1)/2])
      
        #############################  Refinement ##############################
        #if Rinter.shape == (3,3,3):
        #    #print Rinter, Rinter.shape
        #    R_interpolating = RegularGridInterpolator((self.nx, self.ny, self.nz), Rinter,bounds_error=False,fill_value=2.0)     
        #    res = minimize(R_interpolating,np.array([0.01, 0.01, 0.01]), method='nelder-mead',options={'xtol': 1e-2, 'disp': False}) 
        #    #print res.x                                                                               
        #    # the index starts at 0, so if array length is N, which is from 0 to N-1, and the middle index is (N-1+0)/2
        #    # 'corre_center' is a reference point for the displacement
        #    #print (R_index[0]-corre_center[0]+res.x[0],R_index[1]-corre_center[1]+res.x[1],R_index[2]-corre_center[2]+res.x[2])
        #    return ((R_index[0]-corre_center[0]+res.x[0],R_index[1]-corre_center[1]+res.x[1],R_index[2]-corre_center[2]+res.x[2]))
        #else:
        return ((R_index[0]-corre_center[0],R_index[1]-corre_center[1],R_index[2]-corre_center[2]))
       
    def boundarypadding(self,subset,npad):
        # pad dimension is according to (0,1,2)
        x = np.pad(subset,pad_width=npad, mode='constant', constant_values=0)
        return x
        
    def autosubsetsize(self):
        pass
        
    def subset_ref(self,coordinate,size,ref_ind):
        z = coordinate[0]
        x = coordinate[1]
        y = coordinate[2]
        if size%2 == 1:
            subset_radius = (size-1)/2
        elif size%2 ==0:
            assert  size%2 is 1, "The input search subset must be odd number!"      

        if ref_ind == 1:
            subsety = self.image1[z-subset_radius:z+subset_radius+1, x-subset_radius:x+subset_radius+1, y-subset_radius:y+subset_radius+1]
        elif ref_ind == 2:
            subsety = self.image2[z-subset_radius:z+subset_radius+1, x-subset_radius:x+subset_radius+1, y-subset_radius:y+subset_radius+1]
        else:
            print "unknow identifier for the image set"
        return subsety
  
    def normalize(self,subset):
        "this method provides the normalized the region for FFTCC, reducing noise"
        uplim = subset.max()
        lowlim = subset.min()
        rangelim = uplim - lowlim
        norm_sub = (subset-lowlim)/rangelim
        return norm_sub
        
    def FFT(self,subset1,subset2):
        from numpy import unravel_index
        # npad is a tuple of (n_before, n_after) for each dimension
        # e.g.: npad = ((0,0), (1,2), (2,1))        
        x =subset2 # signal function
        y =subset1 # response function    
        sizeFFT = ((x.shape[0]+y.shape[0]),(x.shape[1]+y.shape[1]),(x.shape[2]+y.shape[2]))   
        #sizeFFT = (64,64,64)
        # the size of FFT needs to be 2^n? e.g. (32,32,32)
        fpad2 = y**2
        gpad2 = x**2      
        # put a test unit here
        from scipy.fftpack import fftn, ifftn
        g = ((ifftn(fftn(np.ones(sizeFFT),sizeFFT)*fftn(gpad2,sizeFFT),sizeFFT))**0.5)
        f = ((ifftn(fftn(np.ones(sizeFFT),sizeFFT)*fftn(fpad2,sizeFFT),sizeFFT))**0.5)
        numer = (ifftn(fftn(x,sizeFFT)*np.conjugate(fftn(y,sizeFFT))))
        #g = (np.fft.irfftn(np.fft.rfftn(np.ones(sizeFFT),sizeFFT)*np.fft.rfftn(gpad2,sizeFFT),sizeFFT))**0.5
        #f = (np.fft.irfftn(np.fft.rfftn(np.ones(sizeFFT),sizeFFT)*np.fft.rfftn(fpad2,sizeFFT),sizeFFT))**0.5
        #numer = np.fft.irfftn(np.fft.rfftn(x,sizeFFT)*np.conjugate(np.fft.rfftn(y,sizeFFT)))
        R_xy0 = np.real(numer/g/f)
        #print R_xy0
        # cross correlation matrix
        R_xy1 = R_xy0.take(range(0,x.shape[0]), axis=0)    
        R_xy2 = R_xy1.take(range(0,x.shape[1]), axis=1)  
        R_xy3 = R_xy2.take(range(0,x.shape[2]), axis=2)
        # extract the region of interest for the cross correlation matrix
        R_index = unravel_index(R_xy3.argmax(),R_xy3.shape) 
        print  R_xy3.max()
        # finding the index of the max correlation matrix
        corre_center = ((x.shape[0]+1)/2,(x.shape[1]+1)/2,(x.shape[2]+1)/2)
        # 'corre_center' is a reference point for the displacement
        return (R_index[0]-corre_center[0],R_index[1]-corre_center[1],R_index[2]-corre_center[2])
        # return the displacement
    
    def add(self,id, val1,val2):
        self._oridic[id] = val1
        self._defdic[id] = val2

    def subspre(self):
        zdim = self.point_xyz.shape[0] 
        print "the number of points is" + str(zdim)   
        for i in range(zdim):
            z = int(self.point_xyz[i][0])
            x = int(self.point_xyz[i][1])
            y = int(self.point_xyz[i][2])
            coordinate = np.array([z+2*self.subsize_search,x+2*self.subsize_search,y+2*self.subsize_search])
            # the axis of coordinator is consistant with the physical planes
            subset1 = self.subset_ref(coordinate,self.subsize_fixed,1)
            #print subset1.max()
            # the reference subset1 sets for the subset window
            # the width is the size of target searching area
            subset2 = self.subset_ref(coordinate,self.subsize_search,2)
            self.add(i,subset1,subset2)
        print "data preparation ready!"
        
    def loopcal(self):
        #from math import floor
        zdim = self.point_xyz.shape[0]    
        #self.data1 = []
        #self.data2 = []
        self.subspre()
        nx = np.linspace(-1,1,3)
        ny = np.linspace(-1,1,3)
        nz = np.linspace(-1,1,3)
        self.grid = np.meshgrid(nx, ny, nz)
        for i in range(zdim):
            self.progress = float(i)/zdim*100 
            #obj.y = obj.y + 1
            print round(self.progress)
            # return the progress percentage of DVC process
            #z = int(self.point_xyz[i][0])
            #x = int(self.point_xyz[i][1])
            #y = int(self.point_xyz[i][2])
            #coordinate = np.array([z+2*self.subsize_search,x+2*self.subsize_search,y+2*self.subsize_search])
            ## the axis of coordinator is consistant with the physical planes
            #subset1 = self.subset_ref(coordinate,self.subsize_fixed,1)
            ##print subset1.max()
            ## the reference subset1 sets for the subset window
            ## the width is the size of target searching area
            #subset2 = self.subset_ref(coordinate,self.subsize_search,2)
            ##print subset1.max()
            ## the target subset2 is the region of interest to find the best correlation     
            ##########################  FFT CC #################################
            #subset1 = self._oridic[i]
            #subset2 = self._defdic[i] 
            #nsubset1 = self.normalize(subset1)
            #nsubset2 = self.normalize(subset2)
            # self.FFT(nsubset1,nsubset2) 
            # normalized regions for FFT is required, but for DCC is not necessary
            #result = self.FFT(subset1,subset2)
            #print result
            #self.correlationz[z][x][y] = result[0] 
            #self.correlationx[z][x][y] = result[1]
            #self.correlationy[z][x][y] = result[2]
            #self.point_z_dis[i] = result[0]
            #self.point_x_dis[i] = result[1]
            #self.point_y_dis[i] = result[2]   
            # calculating the correlated displacement in the region
            ####################################################################
            s1 = self._oridic[i]
            s2 = self._defdic[i]    
            ##########################  DCC ####################################
            result  = self.DCC(s1,s2) 
            
            z = int(self.point_xyz[i][0])
            x = int(self.point_xyz[i][1])
            y = int(self.point_xyz[i][2])
            
            coordinate1 = np.array([z+2*self.subsize_search,
                                    x+2*self.subsize_search,
                                    y+2*self.subsize_search])
                                    
            coordinate2 = np.array([z+result[0]+2*self.subsize_search,
                                    x+result[1]+2*self.subsize_search,
                                    y+result[2]+2*self.subsize_search])                       
            
            # the axis of coordinator is consistant with the physical planes                             
 
            #if i==10:
            
            optiSubset1 = self.subset_ref(coordinate1,self.subsize_fixed+14,1) 
            optiSubset2 = self.subset_ref(coordinate2,self.subsize_fixed+14,2)
            
                #print optiSubset1.shape, optiSubset2.shape, "coordinates"
                #print optiSubset1.mean(), optiSubset2.mean()

                
            result2 = self.klt_sub(coordinate1,coordinate2)[1]
            
            #print result2
           
            ####################################################################
            #if i==20:
            #    multiplier = 3       
            #    x1, x2  = coordinate1[0]-multiplier*(self.subsize_fixed-1)/2, coordinate1[0]+multiplier*(self.subsize_fixed-1)/2+1
            #    y1, y2  = coordinate1[1]-multiplier*(self.subsize_fixed-1)/2, coordinate1[1]+multiplier*(self.subsize_fixed-1)/2+1
            #    z1, z2  = coordinate1[2]-multiplier*(self.subsize_fixed-1)/2, coordinate1[2]+multiplier*(self.subsize_fixed-1)/2+1
            #    subset1 = self.image1[x1:x2, y1:y2, z1:z2]
                
            #    xp1, xp2 = coordinate2[0]-multiplier*(self.subsize_fixed-1)/2, coordinate2[0]+multiplier*(self.subsize_fixed-1)/2+1
            #    yp1, yp2 = coordinate2[1]-multiplier*(self.subsize_fixed-1)/2, coordinate2[1]+multiplier*(self.subsize_fixed-1)/2+1
            #    zp1, zp2 = coordinate2[2]-multiplier*(self.subsize_fixed-1)/2, coordinate2[2]+multiplier*(self.subsize_fixed-1)/2+1
                
            #    subset2 = self.image2[xp1:xp2, yp1:yp2, zp1:zp2]
                
            #    vl = Visulization.DataVisulization(subset1, subset1.mean())
            #    vl.contour3d()	
            #    vl = Visulization.DataVisulization(subset2, subset2.mean())
            #    vl.contour3d() 
            #    break
           
            self.correlationz[z][x][y] = result[0]+result2[0]
            self.correlationx[z][x][y] = result[1]+result2[1]
            self.correlationy[z][x][y] = result[2]+result2[2]
            self.point_z_dis[i] = result[0]+result2[0]
            self.point_x_dis[i] = result[1]+result2[1]
            self.point_y_dis[i] = result[2]+result2[2]
                  
    def interpolation(self,xi):
        #from scipy import interpolate
        inter_z = interpolate.LinearNDInterpolator(self.point_xyz, self.point_z_dis, fill_value=0.0, rescale=True)
        inter_x = interpolate.LinearNDInterpolator(self.point_xyz, self.point_x_dis, fill_value=0.0, rescale=True)
        inter_y = interpolate.LinearNDInterpolator(self.point_xyz, self.point_y_dis, fill_value=0.0, rescale=True)
        print "imported!3"
        self.extrapo_z = inter_z.__call__(xi)
        self.extrapo_x = inter_x.__call__(xi)
        self.extrapo_y = inter_y.__call__(xi)

    def klt_sub(self, image1_coordinator, image2_coordinator):
        """
        sub-voxel accuracy
        This module implements KLT on a 3D images
        """
        def warp_a(p, x1, x2, y1, y2, z1, z2):
            """
            WARP_A - Affine warp the image
            WIMG = WARP_A(IMG, P, DST)
            Warp image IMG to WIMG. DST are the destination points, i.e. the corners
            of the template image. P are the affine warp parameters that project
            DST into IMG.
            P = [p1, p4, p7, p10
                p2, p5, p8, p11
                p3, p6, p9, p12];

            Convert affine warp parameters into 3 x 3 warp matrix
            NB affine parameterised as [1 + p1,    p4,     p7,     p10; 
                                    p2,     1 + p5,   p8,     p11;
                                    p3,       p6,     p9+1,   p12]
            """
            Mt = np.concatenate((p, np.array([[0,0,0,1]])), axis=0)
            Mt[0,0] = Mt[0,0] + 1
            Mt[1,1] = Mt[1,1] + 1 
            Mt[2,2] = Mt[2,2] + 1        
            mx, my, mz = np.meshgrid(np.linspace(x1,x2, x2-x1+1), np.linspace(y1,y2, y2-y1+1), np.linspace(z1, z2, z2-z1+1))
            # meshgrid change the orientation of x,y: x becomes the horizontal axis, y becomes the vertical
            # this change would affect the return statement of reshape.
            Cord_xy = np.vstack([mx.flatten(), my.flatten(), mz.flatten()]) # column major vector
            Cord_xy[0,:] = Cord_xy[0,:] 
            Cord_xy[1,:] = Cord_xy[1,:] 
            Cord_xy[2,:] = Cord_xy[2,:] 
            xy = np.concatenate((Cord_xy, np.ones([1,Cord_xy.shape[1]])), axis=0)
            uv = np.dot(Mt,xy)
            # Remove homogeneous
            uv = uv[0:3,:].T
            intermid =  self.fff(uv).reshape(((y2-y1+1),(x2-x1+1),(z2-z1+1)))
            return np.swapaxes(intermid,0,1) # swap the axis to align the image
    
        def jacobian_a(nz, ny, nx):
            # JACOBIAN_A - Compute Jacobian for affine warp
            # DW_DP = JACOBIAN_A(WIDTH, HEIGHT)
            """
            P = [p1, p4, p7, p10
                p2, p5, p8, p11
                p3, p6, p9, p12];
    
            Convert affine warp parameters into 3 x 3 warp matrix
            NB affine parameterised as [1 + p1,    p4,     p7,     p10; [x]
                                        p2,     1 + p5,   p8,     p11; [y]
                                        p3,       p6,     p9+1,   p12  [z]
                                        0         0        0       1 ] [1]
                                                                    
                                    [x,0,0,y,0,0,z,0,0,1,0,0]
                                    [0,x,0,0,y,0,0,z,0,0,1,0]
                                    [0,0,x,0,0,y,0,0,z,0,0,1]                                
            """
            jac_x = np.zeros([nz,ny,nx])
            jac_y = np.zeros([nz,ny,nx])
            jac_z = np.zeros([nz,ny,nx])
            for i in range(nz):
                for j in range(ny):
                    for k in range(nx):
                        jac_x[i,j,k] = k
                        jac_y[i,j,k] = j
                        jac_z[i,j,k] = i
              
            Np = 12
            jac_zero = np.zeros([nz, ny, nx])
            jac_one = np.ones([nz, ny, nx])    
            dW_dp = np.zeros([nz,3*ny,Np*nx]) # 3 is for three dimensions, each dimensional jacobian
            #####################  dW_dp1  ###########################
            dW_dp[0:nz,0:ny,0:nx] = jac_x 
            dW_dp[0:nz,0:ny,0+nx:nx+nx] = jac_zero
            dW_dp[0:nz,0:ny,0+2*nx:nx+2*nx] = jac_zero
            dW_dp[0:nz,0:ny,0+3*nx:nx+3*nx] = jac_y
            dW_dp[0:nz,0:ny,0+4*nx:nx+4*nx] = jac_zero
            dW_dp[0:nz,0:ny,0+5*nx:nx+5*nx] = jac_zero
            dW_dp[0:nz,0:ny,0+6*nx:nx+6*nx] = jac_z
            dW_dp[0:nz,0:ny,0+7*nx:nx+7*nx] = jac_zero
            dW_dp[0:nz,0:ny,0+8*nx:nx+8*nx] = jac_zero
            dW_dp[0:nz,0:ny,0+9*nx:nx+9*nx] = jac_one
            dW_dp[0:nz,0:ny,0+10*nx:nx+10*nx] = jac_zero
            dW_dp[0:nz,0:ny,0+11*nx:nx+11*nx] = jac_zero
            #####################  dW_dp2  ###########################
            dW_dp[0:nz,0+ny:ny+ny,0:nx] = jac_zero
            dW_dp[0:nz,0+ny:ny+ny,0+nx:nx+nx] = jac_x
            dW_dp[0:nz,0+ny:ny+ny,0+2*nx:nx+2*nx] = jac_zero
            dW_dp[0:nz,0+ny:ny+ny,0+3*nx:nx+3*nx] = jac_zero
            dW_dp[0:nz,0+ny:ny+ny,0+4*nx:nx+4*nx] = jac_y
            dW_dp[0:nz,0+ny:ny+ny,0+5*nx:nx+5*nx] = jac_zero
            dW_dp[0:nz,0+ny:ny+ny,0+6*nx:nx+6*nx] = jac_zero
            dW_dp[0:nz,0+ny:ny+ny,0+7*nx:nx+7*nx] = jac_z
            dW_dp[0:nz,0+ny:ny+ny,0+8*nx:nx+8*nx] = jac_zero
            dW_dp[0:nz,0+ny:ny+ny,0+9*nx:nx+9*nx] = jac_zero
            dW_dp[0:nz,0+ny:ny+ny,0+10*nx:nx+10*nx] = jac_one
            dW_dp[0:nz,0+ny:ny+ny,0+11*nx:nx+11*nx] = jac_zero
            ####################  dW_dp3  ###########################
            dW_dp[0:nz,0+2*ny:ny+2*ny,0:nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+nx:nx+nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+2*nx:nx+2*nx] = jac_x
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+3*nx:nx+3*nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+4*nx:nx+4*nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+5*nx:nx+5*nx] = jac_y
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+6*nx:nx+6*nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+7*nx:nx+7*nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+8*nx:nx+8*nx] = jac_z
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+9*nx:nx+9*nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+10*nx:nx+10*nx] = jac_zero
            dW_dp[0:nz,0+2*ny:ny+2*ny,0+11*nx:nx+11*nx] = jac_one  
            return dW_dp
    
        def sd_images(dW_dp, nabla_Ix, nabla_Iy, nabla_Iz, N_p, nz, ny, nx):
            VI_dW_dp = np.zeros([nz, ny, N_p*nx])
            #print dW_dp.shape, nz,ny,nx
            for p in range(N_p):
       	        Tx = nabla_Ix * dW_dp[ 0:nz, 0:ny,         (p*nx)+0:(p*nx)+nx ]
       	        Ty = nabla_Iy * dW_dp[ 0:nz, ny:2*ny,      (p*nx)+0:(p*nx)+nx ]
       	        Tz = nabla_Iz * dW_dp[ 0:nz, 2*ny:3*ny,    (p*nx)+0:(p*nx)+nx ]
       	        VI_dW_dp[:,:,(p*nx)+0:(p*nx)+nx] = (Tx + Ty + Tz)
            return VI_dW_dp
        
        def hessian(VI_dW_dp, N_p, nx):
            # HESSIAN - Compute Hessian
            # H = HESSIAN(VI_DW_DP, N_P, W)
            H = np.zeros([N_p, N_p])
            for i in range(N_p):
   	        h1 = VI_dW_dp[:,:,(i*nx)+0:(i*nx)+nx]
   	        for j in range(N_p):
   	            h2 = VI_dW_dp[:,:,(j*nx)+0:(j*nx)+nx]
                    H[i,j] = np.sum(h1*h2)
   
            #plt.figure(1),imshow(inv(H),cmap='Greys_r'),show()
            return H
    
        def sd_update(VI_dW_dp, error_img, N_p, nx):
            sd_delta_p = np.zeros([N_p,  1])
            for p in range(N_p):
   	        h1 = VI_dW_dp[:,:,(p*nx)+0:(p*nx)+nx]
   	        sd_delta_p[p] = np.sum(h1*error_img)
            return sd_delta_p
   	
        def update_step(warp_p, delta_p):
            # Compute and apply the update
            delta_p = np.reshape(delta_p, (3, 4))
            # Convert affine notation into usual Matrix form - NB transposed
            delta_M = np.concatenate((delta_p, np.array([[0,0,0,1]])), axis=0)
            delta_M[0,0] = delta_M[0,0] + 1
            delta_M[1,1] = delta_M[1,1] + 1
            delta_M[2,2] = delta_M[2,2] + 1
            # Invert compositional warp
            delta_M = inv(delta_M)    
            # Current warp
            warp_M = np.concatenate((warp_p, np.array([[0,0,0,1]])), axis=0)   
            warp_M[0,0] = warp_M[0,0] + 1
            warp_M[1,1] = warp_M[1,1] + 1
            warp_M[2,2] = warp_M[2,2] + 1
            # Compose
            comp_M = np.dot(warp_M, delta_M)
            warp_p = comp_M[0:3,:]
            warp_p[0,0] = warp_p[0,0] - 1
            warp_p[1,1] = warp_p[1,1] - 1
            warp_p[2,2] = warp_p[2,2] - 1	
            return warp_p    
                
        multiplier = 1  
        x1, x2  = image1_coordinator[0]-multiplier*(self.subsize_fixed-1)/2, image1_coordinator[0]+multiplier*(self.subsize_fixed-1)/2+1
        y1, y2  = image1_coordinator[1]-multiplier*(self.subsize_fixed-1)/2, image1_coordinator[1]+multiplier*(self.subsize_fixed-1)/2+1
        z1, z2  = image1_coordinator[2]-multiplier*(self.subsize_fixed-1)/2, image1_coordinator[2]+multiplier*(self.subsize_fixed-1)/2+1
                                    
        xp1, xp2 = image2_coordinator[0]-multiplier*(self.subsize_fixed-1)/2, image2_coordinator[0]+multiplier*(self.subsize_fixed-1)/2+1
        yp1, yp2 = image2_coordinator[1]-multiplier*(self.subsize_fixed-1)/2, image2_coordinator[1]+multiplier*(self.subsize_fixed-1)/2+1
        zp1, zp2 = image2_coordinator[2]-multiplier*(self.subsize_fixed-1)/2, image2_coordinator[2]+multiplier*(self.subsize_fixed-1)/2+1
                
                         
        #image2 = warp_a(p,x1, x2-1, y1, y2-1, z1, z2-1)
        #print x1,x2,y1,y2,z1,z2
    
        points = np.array([[x1,x1,x1,x1,x2,x2,x2,x2],[y1,y1,y2,y2,y1,y1,y2,y2],[z1,z2,z1,z2,z1,z2,z1,z2]])
        points = np.concatenate((points, np.ones([1,points.shape[1]])), axis=0)

        N_p = 12
        n_iters = 15
        
        img = self.image1[x1:x2, y1:y2, z1:z2]
        tmplt = self.image2[xp1:xp2, yp1:yp2, zp1:zp2]
        
        #print img.shape
        #print tmplt.shape
    
        
        #print img.shape,tmplt.shape
    
        nabla_Tx = ndimage.sobel(tmplt,0)
        nabla_Ty = ndimage.sobel(tmplt,1)
        nabla_Tz = ndimage.sobel(tmplt,2)
    
        nz,ny,nx = tmplt.shape
        # 4) Evaluate Jacobian - constant for affine warps
        dW_dp = jacobian_a(nz, ny, nx)
        # 5) Compute steepest descent images, VT_dW_dp
        VT_dW_dp = sd_images(dW_dp, nabla_Tx, nabla_Ty, nabla_Tz, N_p, nz, ny, nx)
        # 6) Compute Hessian and inverse
        H = hessian(VT_dW_dp, N_p, nx)
        #H_gn = hessian(VT_dW_dp, N_p, w);
        #H = np.diag(np.diag(H_gn), 0);
        #H = np.diag(np.diag(H), 0)
        d2G_dp2 = H
        #H_inv = inv(H)
        # Baker-Matthews, Inverse Compositional Algorithm -------------------------
        # NB affine parameterised as [1 + p1, p3, p5; p2, 1 + p4, p6]
        warp_p = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        old_delta_p = np.zeros([N_p,1])
        
        error = np.zeros(n_iters)
        for f in range(n_iters):
   	# 1) Compute warped image with current parameters
   	# Transform into source
             M = np.concatenate((warp_p, np.array([[0,0,0,1]])), axis=0)
             M[0,0] = M[0,0] + 1
             M[1,1] = M[1,1] + 1
             M[2,2] = M[2,2] + 1
             #tpoints = np.dot(M,points).T
       	     IWxp = warp_a(warp_p, x1, x2-1, y1, y2-1, z1, z2-1)
             #print IWxp.shape
             error_img = IWxp - tmplt;
             sd_delta_p = sd_update(VT_dW_dp, error_img, N_p, nx)
           	
             dG_dp = sd_update(VT_dW_dp, error_img, N_p, nx).T
             c = np.dot(dG_dp, dG_dp.T)/np.dot(np.dot(dG_dp, d2G_dp2), dG_dp.T)
           	
             delta_p = c * dG_dp.T
             p_error = delta_p - old_delta_p
             old_delta_p = delta_p
            
             warp_p = update_step(warp_p, delta_p)
             #print np.sqrt(np.sum(p_error**2)) # the newer
             #print np.sqrt(np.sum(error_img**2))  
             error[f] = np.sqrt(np.sum(error_img**2))
                     
             print np.sqrt(np.sum(error_img**2))   
             if(error[f]-error[f-1])>0 and f>=1:
                break
                print "none converged trend!!!" 	 
         
        #vl = Visulization.DataVisulization(IWxp, IWxp.mean())
        #vl.contour3d()	
        #vl = Visulization.DataVisulization(tmplt, tmplt.mean())
        #vl.contour3d() 
        
        	
        
        tpoints = np.dot(M,points).T
        translation = np.mean(tpoints.T,1)-np.mean(points,1)
        print translation
        return warp_p[:,2],translation    
        
if __name__ == '__main__':        
    # need to build a test unit 
    print "This is a test" 
    #import CorriFuc 
    import time
    import scipy.io as sio
    #from LoadDicom import ReadDicom
    #from DataVisu import DataVisulization
    #from PoissonDisk import DiskSamples
    
    import InputOutput 
    import SpatialSampling
    import Visulization
    import CorrelationFunction
    
    from mayavi import mlab
    import mayavi
    
    mat_contents1 = sio.loadmat('/Users/junchaowei/Desktop/Python_DVC2/UP_Research/793-crop-LC-smoothed.mat')
    mat_contents2 = sio.loadmat('/Users/junchaowei/Desktop/Python_DVC2/UP_Research/808-crop-LC-smoothed.mat')

    image1 = mat_contents1['volume']
    image2 = mat_contents2['volume']

    #dicom_volume1 = ReadDicom(PathDicom1)
    #dicom_volume2 = ReadDicom(PathDicom2)
    #dicom_volume1.loadfiles()
    #dicom_volume2.loadfiles()
    #image1 = dicom_volume1.DicArray
    #image2 = dicom_volume2.DicArray
    
    newimage = np.zeros_like(image1)
    threshold = 8000 # choosing the separated level
    
    newimage[image1 < threshold ] = 0
    newimage[image1 >= threshold] = 1
    
    ndimage.binary_opening(newimage, structure=np.ones((3,3,3))).astype(np.int)
    # the following line consumes a lot of computation time@
    newimage2_dist = ndimage.distance_transform_edt(newimage)

    Samples = SpatialSampling.DiskSamples(image1, image2, image1.mean())
    Samples.GenSamples(40,20000)
    point_arr = Samples.samples
    
    import matplotlib.pyplot as plt  
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot(new_point_arr2T[0,:], new_point_arr2T[1,:], new_point_arr2T[2,:],
    #    'o', markersize=10, color='yellow', alpha=0.1)      
    ax.plot(point_arr[0,:], point_arr[1,:], point_arr[2,:],
       'o', markersize=8, color='green', alpha=0.2)

    tpoint_arr = np.transpose(point_arr)
   
    t = time.time()
    ai = Correlation(tpoint_arr,image1,image2,15,23)
    ai.loopcal()  # object result
    elapsed = time.time() - t
    print(elapsed)

    Visulization.DataVisulization(image1,4000).scatterplot(point_arr)
    tpoint_arr = np.transpose(point_arr)
    
    disp = ai.point_x_dis
    
    mesh2 = SpatialSampling.MeshDelaunay(point_arr)
    mesh2.mesh()
    mesh2.alpha_shape(40)
    mesh2.mesh_triangle()
        
    for dummy_i in range(1):   
        Mesh_optimizer = SpatialSampling.MeshOptimizer(point_arr,mesh2.ntri)
        Mesh_optimizer.edge_connector_smooth()
            
    ################  smooth the neigbour displacement point using the neighbor average #########
    edge1 = np.sort(np.array([mesh2.ntri[:,0],mesh2.ntri[:,1]]),axis = 0)
    edge2 = np.sort(np.array([mesh2.ntri[:,0],mesh2.ntri[:,2]]),axis = 0)
    edge3 = np.sort(np.array([mesh2.ntri[:,0],mesh2.ntri[:,3]]),axis = 0)
    edge4 = np.sort(np.array([mesh2.ntri[:,1],mesh2.ntri[:,2]]),axis = 0)
    edge5 = np.sort(np.array([mesh2.ntri[:,1],mesh2.ntri[:,3]]),axis = 0)
    edge6 = np.sort(np.array([mesh2.ntri[:,2],mesh2.ntri[:,3]]),axis = 0)
    edge = np.concatenate((edge1,edge2,edge3,edge4,edge5,edge6),axis=1)
    trans_edge = edge.T
    b = np.ascontiguousarray(trans_edge).view(np.dtype((np.void, trans_edge.dtype.itemsize * trans_edge.shape[1])))
    _, idx = np.unique(b, return_index=True)
    mesh2.edge = trans_edge[idx] # this is the unique edges
    num_smoothes = 0
    for nu in range(num_smoothes):
        for i in range(mesh2.points.shape[0]):
            connector1 = mesh2.edge[mesh2.edge[:,0] == i,1] # find the connectivity from the first column of edge
            connector2 = mesh2.edge[mesh2.edge[:,1] == i,0] # find the connectivity from the second column of edge
            if connector1.size !=0  and connector2.size != 0:
                connector = np.concatenate((connector1,connector2))
            elif connector1.size == 0:
                connector = connector2
            elif connector2.size == 0:
                connector = connector1
            connecting_points = point_arr[connector,:]
            connecting_xdis = disp[connector]
            #print connecting_xdis.shape
            #connecting_ydis = extrapo_y[connector,:]
            #connecting_zdis = extrapo_z[connector,:]
            if connector.shape[0] != 0:
            ##############  uniform smoothing #############
                disp[i] = (np.sum(connecting_xdis)+disp[i])/(connector.shape[0]+1) # plus 'one' is the central displace point 
                #disp[i] = np.median(connecting_xdis)
    # the colors should have the same size and shape as array: "Mesh_optimizer.points.T"
    colors = disp
    @mlab.show
    def main():
        mesh2.view(mesh2.polydata(mesh2.points.T,colors))  
    main()
    mayavi.mlab.colorbar(object=None, title="U(X)(voxel)", orientation="vertical", nb_labels=10, nb_colors=None, label_fmt=None)
    "colorbar reference: http://docs.enthought.com/mayavi/mayavi/auto/mlab_decorations.html"

    plot = Visulization.DataVisulization(image1,image1.mean())
    plot.contour3d()
    plot = Visulization.DataVisulization(image2,image2.mean())
    plot.contour3d()
    surf = plot.surf_points
    print surf.shape
    ai.interpolation(surf)
    plot.plot3dpoint(ai.extrapo_y)


   

            
            
                
                    



    
    
    
    

        
        
        
        
        
    
    
