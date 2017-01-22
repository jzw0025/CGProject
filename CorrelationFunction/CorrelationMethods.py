"""
this module rewrites the Correlation module, which only has a method.

"""

import numpy as np
import DccFort as FDCC
from numpy import unravel_index
from numpy.linalg import inv
from scipy.interpolate import RegularGridInterpolator
import scipy.ndimage as ndimage
from scipy import linalg
from scipy import spatial
import scipy
from numba import jit
from skimage.feature import match_template

def DCC(subset1,subset2,psubset2):
    
    Z, X, Y = subset1.shape
    Zf, Xf, Yf = subset2.shape
    # get the dimension of two subsets
    Z2 = Z // 2
    X2 = X // 2
    Y2 = Y // 2
    npad = ((Z2,Z2),(X2,X2),(Y2,Y2))
    #psubset2 = np.pad(subset2,pad_width=npad, mode='constant', constant_values=0)
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

    R_index = unravel_index(forrij.argmin(),forrij.shape) 
             
    corre_center = ([(Zf-1)/2,(Xf-1)/2,(Yf-1)/2])
      
    return ((R_index[0]-corre_center[0],R_index[1]-corre_center[1],R_index[2]-corre_center[2])), forrij.min() # return the minimum correlation value

def boundarypadding(subset,npad):
        # pad dimension is according to (0,1,2)
        x = np.pad(subset,pad_width=npad, mode='constant', constant_values=0)
        return x 

def calculateSubset(seg_arr1, seg_arr2, point_xyz, subsize, subsearch):
    """
    this method calculates corresponding subsets for given points
    subsize --- the size of whole window size, instead of radius
    subsearch --- the radius whole window size, instead of radius
    """
    ########## help methods ###########   
    
    subsize_fixed = subsize
    subsize_search = subsearch
    
    image1 = boundarypadding(seg_arr1,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
    image2 = boundarypadding(seg_arr2,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
    
    def subset_ref(coordinate,size,ref_ind):
        z = coordinate[0]
        x = coordinate[1]
        y = coordinate[2]
        if size%2 == 1:
            subset_radius = (size-1)/2
        elif size%2 ==0:
            assert  size%2 is 1, "The input search subset must be odd number!"      

        if ref_ind == 1:
            subsety = image1[z-subset_radius:z+subset_radius+1, x-subset_radius:x+subset_radius+1, y-subset_radius:y+subset_radius+1]
        elif ref_ind == 2:
            subsety = image2[z-subset_radius:z+subset_radius+1, x-subset_radius:x+subset_radius+1, y-subset_radius:y+subset_radius+1]
        else:
            print "unknow identifier for the image set"
        return subsety
       
    zdim = point_xyz.shape[0] 
    print "the number of points is" + str(zdim) 
    oridic = {}
    defdic = {}  
    pdefdic = {}
    correlation = {}
    displacement = {}
    for i in range(zdim):
        z = int(point_xyz[i][0])
        x = int(point_xyz[i][1])
        y = int(point_xyz[i][2])
        coordinate = np.array([z+2*subsize_search,x+2*subsize_search,y+2*subsize_search])
        # the axis of coordinator is consistant with the physical planes
        subset1 = subset_ref(coordinate, subsize_fixed, 1)
        #print subset1.max()
        # the reference subset1 sets for the subset window
        # the width is the size of target searching area
        subset2 = subset_ref(coordinate, subsize_search, 2) 
        eachpad = (subsize_fixed//2) # each direction has a pad "(subsize_fixed//2)"
        psubset2 =  subset_ref(coordinate, subsize_search+2*eachpad, 2)  # (ubsize_search+2*eachpad) must be an odd number!
        
        oridic[i] = subset1
        defdic[i] = subset2
        pdefdic[i] = psubset2
        
        result = match_template(psubset2,subset1)
        ij = np.unravel_index(np.argmax(result), result.shape)
        centerx, centery, centerz = (result.shape[0]-1)/2, (result.shape[1]-1)/2, (result.shape[2]-1)/2
        displacement[i] = (ij[0]-centerx, ij[1]-centery, ij[2]-centerz)
        correlation[i] = 1.0 
        #displacement[i], correlation[i] = DCC(subset1, subset2, psubset2)
        print "The Processed Kernal: " + str(i*1.0/zdim*100) #"correlation value is: " + str(correlation[i])
        
    print "data preparation ready!"
    
    return oridic, defdic, correlation, displacement, pdefdic

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
    
              
def KLT(seg_arr1, seg_arr2, image1_coordinator, image2_coordinator, subsize_fixed, iniWarp):

    """
    sub-voxel accuracy
    This module implements KLT on a 3D images
    
    """
    subsearch = 61 # need to be chosen as a number neither too big nor too small
    image1_coordinator = image1_coordinator + 2*subsearch
    image2_coordinator = image2_coordinator + 2*subsearch
    
    image1 = boundarypadding(seg_arr1,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
    image2 = boundarypadding(seg_arr2,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))

    x = np.arange(0, image1.shape[0])
    y = np.arange(0, image1.shape[1])
    z = np.arange(0, image1.shape[2])
        
    fff = RegularGridInterpolator((x, y, z), image1) # interpolating the original image1

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
        
        if (x2-x1+1 != subsize_fixed) or (y2-y1+1 != subsize_fixed) or (z2-z1+1 != subsize_fixed):
            print x2-x1+1, y2-y1+1, z2-z1+1 
        # must make sure return pixels should be agree with the template
                
        mx, my, mz = np.meshgrid(np.linspace(x1,x2, x2-x1+1), np.linspace(y1,y2, y2-y1+1), np.linspace(z1, z2, z2-z1+1))
        # meshgrid change the orientation of x,y: x becomes the horizontal axis, y becomes the vertical
        # this change would affect the return statement of reshape.
        Cord_xy = np.vstack([mx.flatten(), my.flatten(), mz.flatten()]) # column major vector
        #Cord_xy[0,:] = Cord_xy[0,:] 
        #Cord_xy[1,:] = Cord_xy[1,:] 
        #Cord_xy[2,:] = Cord_xy[2,:] 
        xy = np.concatenate((Cord_xy, np.ones([1,Cord_xy.shape[1]])), axis=0)
        uv = np.dot(Mt,xy)
        # Remove homogeneous
        uv = uv[0:3,:].T
        #if uv.max()>=(y2-y1+1) or uv.max()>=(x2-x1+1) or uv.max()>=(z2-z1+1):
        #    print "the maxier is: " +str(uv.max())
        
        if uv.min()<=0:
            print "the miner is:" + str(uv.min())
        
        intermid =  fff(uv).reshape(((y2-y1+1),(x2-x1+1),(z2-z1+1)))
        return np.swapaxes(intermid,0,1) # swap the axis to align the image
                
    multiplier = 1  
    x1, x2  = image1_coordinator[0]-multiplier*(subsize_fixed-1)/2, image1_coordinator[0]+multiplier*(subsize_fixed-1)/2+1
    y1, y2  = image1_coordinator[1]-multiplier*(subsize_fixed-1)/2, image1_coordinator[1]+multiplier*(subsize_fixed-1)/2+1
    z1, z2  = image1_coordinator[2]-multiplier*(subsize_fixed-1)/2, image1_coordinator[2]+multiplier*(subsize_fixed-1)/2+1
                                    
    xp1, xp2 = image2_coordinator[0]-multiplier*(subsize_fixed-1)/2, image2_coordinator[0]+multiplier*(subsize_fixed-1)/2+1
    yp1, yp2 = image2_coordinator[1]-multiplier*(subsize_fixed-1)/2, image2_coordinator[1]+multiplier*(subsize_fixed-1)/2+1
    zp1, zp2 = image2_coordinator[2]-multiplier*(subsize_fixed-1)/2, image2_coordinator[2]+multiplier*(subsize_fixed-1)/2+1
                                        
    points = np.array([[x1,x1,x1,x1,x2,x2,x2,x2],[y1,y1,y2,y2,y1,y1,y2,y2],[z1,z2,z1,z2,z1,z2,z1,z2]])
    points = np.concatenate((points, np.ones([1,points.shape[1]])), axis=0)
    
    #img = image1[x1:x2, y1:y2, z1:z2]
    tmplt = image2[xp1:xp2, yp1:yp2, zp1:zp2]
    
    nabla_Tx = ndimage.sobel(tmplt,0)
    nabla_Ty = ndimage.sobel(tmplt,1)
    nabla_Tz = ndimage.sobel(tmplt,2)
    
    N_p = 12
    n_iters = 10
    
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
    # warp_p = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    warp_p = iniWarp # assign a initial value
    old_delta_p = np.zeros([N_p,1])
    condition = False
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
    
             ###########	
       	     #Iwxp_mean = np.mean(IWxp)
             #tmplt_mean = np.mean(tmplt)
       	     #delta_Iwxp = np.sqrt((np.sum(IWxp-Iwxp_mean))**2)
             #delta_tmplt = np.sqrt((np.sum(tmplt-tmplt_mean))**2)+np.finfo(float).eps # avoid dividing by zero
             #print "delta_Iwxp is: "+str(delta_tmplt)
             #error_img = (IWxp-Iwxp_mean) - delta_Iwxp/delta_tmplt*(tmplt-tmplt_mean)	
	     #error_img = IWxp - tmplt;
             #sd_delta_p = sd_update(VT_dW_dp, error_img, N_p, nx)
             #dG_dp = sd_update(VT_dW_dp, error_img, N_p, nx).T
             #c = np.dot(dG_dp, dG_dp.T)/np.dot(np.dot(dG_dp, d2G_dp2), dG_dp.T)
             #delta_p = c * dG_dp.T
	     
	     ############
	     #print "new system test!"
	     IWxp_mean = np.mean(IWxp)
	     tmplt_mean = np.mean(tmplt)
	     delta_Iwxp = np.sqrt((np.sum(IWxp-IWxp_mean)**2))+np.finfo(float).eps # avoid dividing by zero
	     delta_tmplt = np.sqrt((np.sum(tmplt-tmplt_mean)**2))+np.finfo(float).eps # avoid dividing by zero
	     error_img = (IWxp-IWxp_mean)/delta_Iwxp - (tmplt-tmplt_mean)/delta_tmplt 
	     dG_dp = sd_update(VT_dW_dp, error_img, N_p, nx).T
             c = np.dot(dG_dp, dG_dp.T)/np.dot(np.dot(dG_dp, d2G_dp2), dG_dp.T)
             
             #delta_p = -c * dG_dp.T
             
             delta_p = c * dG_dp.T*delta_tmplt
             
             p_error = delta_p - old_delta_p
             old_delta_p = delta_p
            
             warp_p = update_step(warp_p, delta_p)
            
             error_img2 = IWxp - tmplt
             
             error[f] = np.sqrt(np.sum(error_img2**2))
                     
             print "the least square erros is: "+str(np.sqrt(np.sum(error_img2**2))) 
              
             if(error[f]-error[f-1])>=0 and f>=1:
                if f < 4:
                    condition = True
                    print "none converged trend!!!" 
                break
                	 
         
        #vl = Visulization.DataVisulization(IWxp, IWxp.mean())
        #vl.contour3d()	
        #vl = Visulization.DataVisulization(tmplt, tmplt.mean())
        #vl.contour3d() 
   
    tpoints = np.dot(M,points).T
    translation = np.mean(tpoints.T,1)-np.mean(points,1)
    print "the translation vector is:"+str(translation)
    return warp_p, translation, IWxp, tmplt, condition
    
def ConnectMesh(ntri, pointNumber):
    """
    this function returns the connectivity of triangle mesh.
    
    """
    edge1 = np.sort(np.array([ntri[:,0],ntri[:,1]]),axis = 0) # making all possible between the edges
    edge2 = np.sort(np.array([ntri[:,0],ntri[:,2]]),axis = 0)
    edge3 = np.sort(np.array([ntri[:,0],ntri[:,3]]),axis = 0)
    edge4 = np.sort(np.array([ntri[:,1],ntri[:,2]]),axis = 0)
    edge5 = np.sort(np.array([ntri[:,1],ntri[:,3]]),axis = 0)
    edge6 = np.sort(np.array([ntri[:,2],ntri[:,3]]),axis = 0)
    
    edge = np.concatenate((edge1,edge2,edge3,edge4,edge5,edge6),axis=1)
    trans_edge = edge.T
    
    b = np.ascontiguousarray(trans_edge).view(np.dtype((np.void, trans_edge.dtype.itemsize * trans_edge.shape[1])))
    _, idx = np.unique(b, return_index=True) # find unique set
    
    edge = trans_edge[idx] # this is the unique edges
    
    connector1 = edge[ edge[:,0] == pointNumber,1 ] # find the connectivity from the first column of edge
    connector2 = edge[ edge[:,1] == pointNumber,0 ] # find the connectivity from the second column of edge
    
    if connector1.size !=0  and connector2.size != 0:
        connector = np.concatenate((connector1,connector2))
        
    elif connector1.size == 0:
        connector = connector2
        
    elif connector2.size == 0:
        connector = connector1
    else:
        print "there is no connecting point to the target!"
        
    return connector # the number of point/node connected to the point
    
def SVD(target, source):
    """
    the input data:
    target --- N-D array e.g. shape: (1000, 3)
    source --- N-D array e.g. shape: (1000, 3) 
    the return data:
            
    """
    targetTree = spatial.KDTree(target) # create K-D tree for the input target
            
    ######### the following code goes iterately ##########
    translation = np.zeros([1,3])
    
    ##########   translation  ########
    vector = target.mean(0) - source.mean(0)
    source = source + vector
    ##########  rotation  ############
    sourcePair_list = []
    targetPair_list = []
                
    for i in range(source.shape[0]):
        # find the closest point in the tree-data
        results = targetTree.query(source[i])
        sourcePair_list.append(list(source[i]))
        targetPair_list.append(list(target[results[1]]))
                    
    matrix1 = np.matrix(sourcePair_list)
    matrix2 = np.matrix(targetPair_list)
    # building the corresponding matrix
    resultMatrix = matrix1.T.dot(matrix2)
    # perform the SVD algorithm
    U, s, Vh = linalg.svd(resultMatrix)
    Rotation = Vh.T.dot(U.T)
    source2 = Rotation.dot(np.matrix(source).T).T
    #print newsource == newsource2
    source = np.array(source2) # updating the source points
    translation = translation + vector
            
    return source, Rotation, translation 
    
def fixed_SVD(target, source):
    """
    this module calculates the fixed SVD for the paired input points
    
    """
    translation = np.zeros([1,3])
    vector = target.mean(0) - source.mean(0) # first dimension mean
    source = source + vector
    ##########  rotation  ############
    matrix1 = np.matrix(source)
    matrix2 = np.matrix(target)
    # building the corresponding matrix
    resultMatrix = matrix1.T.dot(matrix2)
    # perform the SVD algorithm
    U, s, Vh = linalg.svd(resultMatrix)
    Rotation = Vh.T.dot(U.T)
    source2 = Rotation.dot(np.matrix(source).T).T
    #print newsource == newsource2
    source = np.array(source2) # updating the source points
    translation = translation + vector # translation is zero initially, this does not need to add iteratively, but keep the original formula
            
    return source, Rotation, translation 
    
    
def ICP(image2, point1, point2, iterations):
    change_points = point2
    rotation_matrix_list =  []  
    translation_vector_list = [] 
    for i in range(iterations): # this SVD registration is followed by image registration at each step.
        old_change_points = change_points # save a copy of old points
        change_points, R, V = SVD(point1, change_points) # the translation is not used
        print "the sum sqare for the point cloud is: " + str(np.sqrt(sum(sum((change_points-old_change_points)**2))))
        rotation2 = R
        rotation_matrix_list.append(R)
        #translation = self.mean1_vector.T[0] - np.dot(np.dot(self.eig_vec_cov1,self.eig_vec_cov2.T),self.mean2_vector.T[0])
        translation2 =  ((point2.mean(0)).dot(rotation2)-point1.mean(0)).dot(linalg.inv(rotation2)) # recalculating the translation using calculated rotation
        translation_vector_list.append(translation2)
        newimage2 = ndimage.interpolation.affine_transform(image2,rotation2,order=3,offset=translation2,cval=0.0) # rigid affine registration 
           
    return newimage2, change_points, rotation_matrix_list, translation_vector_list            
              
class CornerDetector:
    """
    this class implements HarrisCorner3D for 3D images.
    ref: http://www.cse.psu.edu/~rtc12/CSE486/lecture06.pdf 
    Date: 7/1/16
        
    """
    def __init__(self,image1, windowSize):
        """ 
        the boundary subset needs to be considered
        windowSize is window radius
        """
        self.subsearch = windowSize # need to be chosen as a number neither too big nor too small
        ############## padding the image1 #############
        self.image1 = boundarypadding(image1,((self.subsearch,self.subsearch),(self.subsearch,self.subsearch),(self.subsearch,self.subsearch)))   
        
        self.Ix = ndimage.sobel(self.image1,0) # pre-calculate the derivatives
        self.Iy = ndimage.sobel(self.image1,1)
        self.Iz = ndimage.sobel(self.image1,2)
        
        self.Ixx = self.Ix*self.Ix
        self.Iyy = self.Iy*self.Iy
        self.Izz = self.Iz*self.Iz
        self.Ixy = self.Ix*self.Iy
        self.Ixz = self.Ix*self.Iz
        self.Iyz = self.Iy*self.Iz
        
        # pre-compute gaussian window
        self.f = self.getGaussianWeight(windowSize*2+1, windowSize*2+1, windowSize*2+1, 1) # Sigma=1 adjust the window to the given subset size
        self.harrisValue = []
        
    def getGaussianWeight(self,x,y,z,Size):
            
        """
        *** This is a Help Function***
                
        this function returns the gaussian weights in the subset form (size, size, size)
                
        #################
        ## 3D Gaussian ##     f(x, y z) = A*exp(-((x-x0)^2/a+(y-y0)^2/b+(z-z0)^2/c))
        #################
                
        """
                
        A = 0.1
        x1,x2,y1,y2,z1,z2 = 0,x-1, 0,y-1, 0,z-1
                
        mx, my, mz = np.meshgrid(np.linspace(x1,x2, x2-x1+1), np.linspace(y1,y2, y2-y1+1), np.linspace(z1, z2, z2-z1+1))
        centerx, centery, centerz = (x-1)/2, (y-1)/2, (z-1)/2
                
        # meshgrid change the orientation of x,y: x becomes the horizontal axis, y becomes the vertical
        # this change would affect the return statement of reshape.
        #Cord_xy = np.vstack([mx.flatten(), my.flatten(), mz.flatten()]) # column major vector
        #print mx-centerx
        f = A*scipy.exp(-((mx-centerx)**2/float(Size) + (my-centery)**2/float(Size) + (mz-centerz)**2/float(Size)))
        
        return f

    def harrisCorner3D(self, points):
        """
        This method implements the Harris Corner Detection Algorithm
        image1 dim: 
        points dim: 1 by 3

        """
        points = points + self.subsearch # offset points in padded image.
        
        # computer x and y derivatives of image
        # from here: everything move from 2D to 3D
        
        window_wide = self.subsearch # the size of window radius
    
        i, j, k =  int(points[0]), int(points[1]), int(points[2])
        
        Sx2 = self.Ixx[ i-window_wide:i+window_wide+1,j-window_wide:j+window_wide+1, k-window_wide:k+window_wide+1] # easy sqared window
        Sy2 = self.Iyy[ i-window_wide:i+window_wide+1,j-window_wide:j+window_wide+1, k-window_wide:k+window_wide+1]
        Sz2 = self.Izz[ i-window_wide:i+window_wide+1,j-window_wide:j+window_wide+1, k-window_wide:k+window_wide+1]
        Sxy = self.Ixy[ i-window_wide:i+window_wide+1,j-window_wide:j+window_wide+1, k-window_wide:k+window_wide+1]
        Sxz = self.Ixz[ i-window_wide:i+window_wide+1,j-window_wide:j+window_wide+1, k-window_wide:k+window_wide+1]
        Syz = self.Iyz[ i-window_wide:i+window_wide+1,j-window_wide:j+window_wide+1, k-window_wide:k+window_wide+1]
        
        fSx2 = np.sum(self.f*Sx2)
        fSy2 = np.sum(self.f*Sy2)
        fSz2 = np.sum(self.f*Sz2)
        fSxy = np.sum(self.f*Sxy)
        fSxz = np.sum(self.f*Sxz)
        fSyz = np.sum(self.f*Syz)
        
        Hxyz = np.array([[fSx2, fSxy, fSxz],[fSxy, fSy2, fSyz],[fSxz, fSyz, fSz2]])
        
        Kcont = 0.25 # k is an empirically determined constant; k = 0.04 - 0.06
        
        #print linalg.det(Hxyz), (np.trace(Hxyz))**2
        
        #R = linalg.det(Hxyz) - Kcont*(np.trace(Hxyz))**2 # calculate the R value
        eig = abs(np.linalg.eigvals(Hxyz))
        print "eig is" +str(eig)
        return eig.min()
        
    def harrisPoint(self, point_arr1, threshold):
        """
        this method ultilizes
        """    
        for hi in range(point_arr1.shape[0]):
            self.harrisValue.append(self.harrisCorner3D(point_arr1[hi,:]))

        counts = 0
        HarrisPoints = []
        
        index = 0
        for harris in self.harrisValue:
            print harris
            if harris > threshold or harris < -threshold:
                counts+=1
                HarrisPoints.append(point_arr1[index,:]) # the points in HarrisPoints is without offset
                index+=1   
                
        print "harris feature ratio is: "+ str(float(counts)/len(self.harrisValue))
        
        HarrisArray = np.ceil(np.array(HarrisPoints))
        
        return HarrisArray
        
    def visuHarris(self,harrisSubPoint):

        harrisRadius = self.subsearch 
        
        harrisSubPoint = harrisSubPoint + harrisRadius
    
        #for isub in range(HarrisArray.shape[0]):
        harrisSubPoint = harrisSubPoint.astype(int)
    
        harrisSubset = self.image1[ harrisSubPoint[0]-harrisRadius:harrisSubPoint[0]+harrisRadius+1,
                                    harrisSubPoint[1]-harrisRadius:harrisSubPoint[1]+harrisRadius+1,
                                    harrisSubPoint[2]-harrisRadius:harrisSubPoint[2]+harrisRadius+1]
        #vl = Visulization.DataVisulization(harrisSubset, 1200)
        #vl.contour3d()
        return harrisSubset
        
def calculateSubset_PCL(seg_arr1, seg_arr2, point_xyz, subsize, subsearch):
    """
    this method calculates corresponding subsets for given points
    
    """
    ########## help methods ###########   
    
    subsize_fixed = subsize
    subsize_search = subsearch
    
    image1 = boundarypadding(seg_arr1,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
    image2 = boundarypadding(seg_arr2,((2*subsearch,2*subsearch),(2*subsearch,2*subsearch),(2*subsearch,2*subsearch)))
    
    def subset_ref(coordinate,size,ref_ind):
        z = coordinate[0]
        x = coordinate[1]
        y = coordinate[2]
        if size%2 == 1:
            subset_radius = (size-1)/2
        elif size%2 ==0:
            assert  size%2 is 1, "The input search subset must be odd number!"      

        if ref_ind == 1:
            subsety = image1[z-subset_radius:z+subset_radius+1, x-subset_radius:x+subset_radius+1, y-subset_radius:y+subset_radius+1]
        elif ref_ind == 2:
            subsety = image2[z-subset_radius:z+subset_radius+1, x-subset_radius:x+subset_radius+1, y-subset_radius:y+subset_radius+1]
        else:
            print "unknow identifier for the image set"
        return subsety
       
    zdim = point_xyz.shape[0] 
    print "the number of points is" + str(zdim) 

    correlation = {}
    displacement = {}
    for i in range(zdim):
        z = int(point_xyz[i][0])
        x = int(point_xyz[i][1])
        y = int(point_xyz[i][2])
        coordinate = np.array([z+2*subsize_search,x+2*subsize_search,y+2*subsize_search])
        # the axis of coordinator is consistant with the physical planes
        subset1 = subset_ref(coordinate, subsize_fixed, 1)
        #print subset1.max()
        # the reference subset1 sets for the subset window
        # the width is the size of target searching area
        subset2 = subset_ref(coordinate, subsize_search, 2) 
        eachpad = (subsize_fixed//2) # each direction has a pad "(subsize_fixed//2)"
        psubset2 =  subset_ref(coordinate, subsize_search+2*eachpad, 2)  # (ubsize_search+2*eachpad) must be an odd number!
        displacement[i], correlation[i] = DCC(subset1, subset2, psubset2)
        print i*1.0/zdim*100 #"correlation value is: " + str(correlation[i])
        
    print "data preparation ready!"
    
    return displacement

def extractTriangle(tetras):
    """
    this method creates triangle elements from the tetrahedron element
    """
    if tetras.shape[1]!=4:
        print "invalid format of input data! The input data show be [N,4] dimensions"
    
    pointA = tetras[:,0]
    pointB = tetras[:,1]
    pointC = tetras[:,2]
    pointD = tetras[:,3]
    
    triangle1 = np.zeros([tetras.shape[0],3]) # making temporary numpy array
    triangle2 = np.zeros([tetras.shape[0],3])
    triangle3 = np.zeros([tetras.shape[0],3])
    triangle4 = np.zeros([tetras.shape[0],3])
    
    ####### triangle 1 (pointA, pointB, pointC) #########
    triangle1[:,0] = pointA
    triangle1[:,1] = pointB
    triangle1[:,2] = pointC
    
    ####### triangle 2 (pointA, pointB, pointD) #########
    triangle2[:,0] = pointA
    triangle2[:,1] = pointB
    triangle2[:,2] = pointD
    
    ####### triangle 3 (pointA, pointC, pointD) #########
    triangle3[:,0] = pointA
    triangle3[:,1] = pointC
    triangle3[:,2] = pointD
    
    ####### triangle 4 (pointB, pointC, pointD) #########
    triangle4[:,0] = pointB
    triangle4[:,1] = pointC
    triangle4[:,2] = pointD
    
    triangle1 = np.sort(triangle1)
    triangle2 = np.sort(triangle2)
    triangle3 = np.sort(triangle3)
    triangle4 = np.sort(triangle4)
    
    c12 = np.concatenate((triangle1,triangle2),axis=0) # concatenating triangle arrays (1,2)
    c123 = np.concatenate((c12,triangle3),axis=0) 
    c1234 = np.concatenate((c123,triangle4),axis=0)
    
    unique_list = {} # this is dictionary for storing the triangle element
    for i in range(c1234.shape[0]):
        tuple_c12345 = tuple(c1234[i,:].tolist())
        if tuple_c12345 not in unique_list.keys():
            unique_list[tuple_c12345] = 1
        else:
            unique_list[tuple_c12345] = unique_list.get(tuple_c12345)+1
    
    out = []
    for j in unique_list.keys():
        if unique_list.get(j) == 1:
            out.append(j)

    return out # create a big set with tuple elements 
    
#def SVD(target, source):
#    """
#    the input data:
#    target --- N-D array e.g. shape: (1000, 3)
#    source --- N-D array e.g. shape: (1000, 3) 
#    the return data:
#    """
#    targetTree = spatial.KDTree(target) # create K-D tree for the input target 
#    ######### the following code goes iterately ##########
#    translation = np.zeros([1,3])
#    ##########   translation  ########
#    vector = target.mean(0) - source.mean(0)
#    source = source + vector
#    ##########  rotation  ############
#    sourcePair_list = []
#    targetPair_list = []       
#    for i in range(source.shape[0]):
#        # find the closest point in the tree-data
#        results = targetTree.query(source[i])
#        sourcePair_list.append(list(source[i]))
#        targetPair_list.append(list(target[results[1]]))           
#    matrix1 = np.matrix(sourcePair_list)
#    matrix2 = np.matrix(targetPair_list)
#    # building the corresponding matrix
#    resultMatrix = matrix1.T.dot(matrix2)
#    # perform the SVD algorithm
#    
#    try:
#        U, s, Vh = linalg.svd(resultMatrix)
#    except ValueError: 
#        return source, np.eye(3), np.zeros((1,3))
#        
#    Rotation = Vh.T.dot(U.T)
#    source2 = Rotation.dot(np.matrix(source).T).T
#    #print newsource == newsource2
#    source = np.array(source2) # updating the source points
#    translation = translation + vector
#    return source, Rotation, translation   
    
def ICP_KLT(point1, point2, iterations):        
    change_points = point2
    translation2 = [0.0,0.0,0.0]
    for i in range(iterations): # this SVD registration is followed by image registration at each step.
        old_change_points = change_points
        change_points, R, V = SVD(point1, change_points)
        print "the sum sqare for the point cloud is: " + str(np.sqrt(sum(sum((change_points-old_change_points)**2))))
        rotation2 = R
        #translation = self.mean1_vector.T[0] - np.dot(np.dot(self.eig_vec_cov1,self.eig_vec_cov2.T),self.mean2_vector.T[0])
        #translation2 = ((point2.mean(0)).dot(rotation2)-point1.mean(0)).dot(linalg.inv(rotation2))   
        translation2 = translation2 + V          
    
    return translation2     
    
def subPoint(point_arr, center, radius):  
               
    x1 = center[0] - radius
    x2 = center[0] + radius 
        
    y1 = center[1] - radius
    y2 = center[1] + radius   
        
    z1 = center[2] - radius
    z2 = center[2] + radius
        
    index_x = np.logical_and(point_arr[:,0]>=x1, point_arr[:,0]<=x2)
    index_y = np.logical_and(point_arr[:,1]>=y1, point_arr[:,1]<=y2)
    index_z = np.logical_and(point_arr[:,2]>=z1, point_arr[:,2]<=z2)    
    index = np.logical_and(np.logical_and(index_x, index_y),index_z)
        
    return point_arr[index]

