def createPoints(image1):
    sx, sy, sz = image1.shape
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