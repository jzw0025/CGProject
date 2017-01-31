"""
this module runs the strain and post processing for the 
Digital Volume Correlation
"""
import numpy as np
from scipy import interpolate

def interp2Strain(points, displacements, sample_length, coordinators):
    """
    this function does the interpolation over the neighbor connecting points
    input:
        points --- (N,3)
        displacements --- (3, N)  List--Array
        sample_length --- float
        coordinators --- (1,3)
    output:
        strains --- (1,6): E11,E22,E33,E12,E13,E23
    """
    displacementsX = displacements[0]
    displacementsY = displacements[1]
    displacementsZ = displacements[2]
    
    U = interpolate.LinearNDInterpolator(points, displacementsX, fill_value=0.0)
    V = interpolate.LinearNDInterpolator(points, displacementsY, fill_value=0.0)
    W = interpolate.LinearNDInterpolator(points, displacementsZ, fill_value=0.0)
    
    xcor = coordinators[0]
    ycor = coordinators[1]
    zcor = coordinators[2]
    
    dudx = (U([xcor+sample_length,ycor,zcor])-U([xcor-sample_length,ycor,zcor]))[0]/2/sample_length
    dudy = (U([xcor,ycor+sample_length,zcor])-U([xcor,ycor-sample_length,zcor]))[0]/2/sample_length
    dudz = (U([xcor,ycor,zcor+sample_length])-U([xcor,ycor,zcor-sample_length]))[0]/2/sample_length

    dvdx = (V([xcor+sample_length,ycor,zcor])-V([xcor-sample_length,ycor,zcor]))[0]/2/sample_length
    dvdy = (V([xcor,ycor+sample_length,zcor])-V([xcor,ycor-sample_length,zcor]))[0]/2/sample_length
    dvdz = (V([xcor,ycor,zcor+sample_length])-V([xcor,ycor,zcor-sample_length]))[0]/2/sample_length
    
    dwdx = (W([xcor+sample_length,ycor,zcor])-W([xcor-sample_length,ycor,zcor]))[0]/2/sample_length
    dwdy = (W([xcor,ycor+sample_length,zcor])-W([xcor,ycor-sample_length,zcor]))[0]/2/sample_length
    dwdz = (W([xcor,ycor,zcor+sample_length])-W([xcor,ycor,zcor-sample_length]))[0]/2/sample_length
    
    F_matrix = np.array([[dudx,dudy,dudz], [dvdx,dvdy,dvdz], [dwdx,dwdy,dwdz]]) + np.eye(3)
    
    Strain =  np.dot(F_matrix.T,F_matrix) - np.eye(3)
    
    StrainList = [Strain[0,0], Strain[1,1], Strain[2,2], Strain[0,1], Strain[0,2],Strain[1,2]]
    
    return StrainList
    
def getEdges(ntri):
    edge1 = np.sort(np.array([ntri[:,0], ntri[:,1]]),axis = 0)
    edge2 = np.sort(np.array([ntri[:,0], ntri[:,2]]),axis = 0)
    edge3 = np.sort(np.array([ntri[:,0], ntri[:,3]]),axis = 0)
    edge4 = np.sort(np.array([ntri[:,1], ntri[:,2]]),axis = 0)
    edge5 = np.sort(np.array([ntri[:,1], ntri[:,3]]),axis = 0)
    edge6 = np.sort(np.array([ntri[:,2], ntri[:,3]]),axis = 0)
    edge = np.concatenate((edge1,edge2,edge3,edge4,edge5,edge6),axis=1)
    trans_edge = edge.T
    b = np.ascontiguousarray(trans_edge).view(np.dtype((np.void, trans_edge.dtype.itemsize * trans_edge.shape[1])))
    _, idx = np.unique(b, return_index=True)
    return trans_edge[idx] # this is the unique edges
    
def getConnector(iP, edges):
    connector1 = edges[edges[:,0] == iP,1] # find the connectivity from the first column of edge
    connector2 = edges[edges[:,1] == iP,0] # find the connectivity from the second column of edge
    if connector1.size !=0  and connector2.size != 0:
        connector = np.concatenate((connector1,connector2))
    elif connector1.size == 0:
        connector = connector2
    elif connector2.size == 0:
        connector = connector1
    else:
        print "the point is isolated without any connecting edges!"
    return connector # np array
    
def getConnectorValue(value, connector):
    """
    input: 
        value --- [3, N] np.array
        connector --- [1,N]
    output:
        connectingDisp --- [3,M]
    """
    if type(value)== list:
        value = np.array(value)
        
    connectingDisp = value[:,connector]
    return connectingDisp
    
def calculateStrain(points, ntri, displacements, sample_length, num_smoothes):
    """
    input:
        ntri --- input element (N,3)
        num_smoothes --- int N
        components  ---- List [[np.array],[np.array],...]  
        "containing the components needed to be calculated"  [X, Y, Z]     displacement     
        points --- (N,3)
    output:
        strain --- (6,N) sublist
    """
    mesh_edge = getEdges(ntri) # this is the unique edges
    s11,s22,s33,s12,s13,s23 = [],[],[],[],[],[]
    for i in range(points.shape[0]):
            print "the computation process for the strain(%)" + str(i*100.0/points.shape[0])
            connector = getConnector(i,mesh_edge)
            connecting_points = points[connector,:]
            coordinators = points[i,:]
            connecting_displacements = getConnectorValue(displacements, connector)
            SCs = interp2Strain(connecting_points, connecting_displacements, sample_length, coordinators) 
            ## smoothing for each strain elements.
            s11.append(SCs[0]) 
            s22.append(SCs[1]) 
            s33.append(SCs[2]) 
            s12.append(SCs[3]) 
            s13.append(SCs[4])
            s23.append(SCs[5])
    return [s11, s22, s33, s12, s13, s23]

def smoothStrains(strain, ntri, num_smoothes):
    """
    input:
        strain --- [6,N] list
        connector ---[1,N] array
        num_smoothes --- int
    output:    
        strain --- [6,N] np array
    """
    if type(strain) == list:
        strain = np.array(strain)
    si = strain.shape
    if len(si) == 1:
        strain_number  = si[0]
    else:
        if si[0] > si[1]:
            strain_number = si[0]
        else:
            strain_number = si[1]
    
    mesh_edge = getEdges(ntri) 
    for nu in range(num_smoothes):
        print strain.shape[1]
        for ni in range(strain_number):  # this loop may be simplied used Ax=B matrix assembly
            print "the computation process for the strain smooth(%)" + str(ni*100.0/strain.shape[1])
            connector = getConnector(ni,mesh_edge)
            np.append(connector,ni) # including itself for smoothing
            connecting_strains =  getConnectorValue(strain, connector)
            strain[:,ni] = np.mean(connecting_strains, 1)

if __name__ == "__main__":
    print "this a unit test process!"
    import SaveLoad
    import SpatialSampling
    import Visualization
    import scipy.ndimage as ndimage
    inputDVC = SaveLoad.DVCdata("/Users/junchaowei/Desktop/HighResolutionStrain/")
    point_arr1 = inputDVC.getPoints1()
    Sample_Density = 12
    
    mesh2 = SpatialSampling.MeshDelaunay(point_arr1.T) # use the original point_arr1 
    mesh2.mesh()
    mesh2.alpha_shape(Sample_Density)
    mesh2.mesh_triangle() # find unique triangles
    
    for dummy_i in range(1):   
        Mesh_optimizer = SpatialSampling.MeshOptimizer(point_arr1.T,mesh2.ntri)
        Mesh_optimizer.edge_connector_smooth()
      
    displacements = [inputDVC.getDisplacementX()[0], 
                                inputDVC.getDisplacementY()[0], 
                                inputDVC.getDisplacementZ()[0]]
    
    strain = calculateStrain(point_arr1, mesh2.ntri, displacements, Sample_Density/2, 2)
    strain = np.array(strain)
    a = np.mean(strain,1)
    smoothStrains(strain, mesh2.ntri,4)
    b =  np.mean(strain,1)
    
    image1 = inputDVC.getImage1()
    plot = Visualization.DataVisulization(image1, image1.mean()+1400)
    plot.contour3d()
    surf = plot.surf_points
    inter_z = interpolate.LinearNDInterpolator(mesh2.points, strain[0,:], fill_value=0.0, rescale=True)        
    extrapo_z = inter_z.__call__(surf)
    plot.plot3dpoint(extrapo_z)
    
#### strain plot ####

    ix, iy, iz = image1.shape
    new_image1 = np.empty((ix+20,iy+20,iz+20))  
    new_image1[20:-20,20:-20,20:-20] = image1[10:-10,10:-10,10:-10]
    
    trans_point = np.empty_like(mesh2.points)
    trans_point[:,0] = mesh2.points[:,0] + 10  
    trans_point[:,1] = mesh2.points[:,1] + 10
    trans_point[:,2] = mesh2.points[:,2] + 10 
    
    plot = Visualization.DataVisulization(new_image1, 8000)
    plot.contour3d()
    surf = plot.surf_points
    inter_z = interpolate.LinearNDInterpolator(trans_point, strain[5,:], fill_value=0.0, rescale=True)        
    extrapo_z = inter_z.__call__(surf)
    plot.plot3dpoint(extrapo_z)
    
    import mayavi
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    