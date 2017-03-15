"""
feature extacted displacements
"""
import numpy as np
from math import *
from mayavi import mlab
from mayavi.mlab import *
import mayavi
from sklearn.neighbors import KDTree
import SpatialSampling
from scipy import interpolate

def minFeaturePairing(d1,d2):
    """
    
    This function matching the features from d1 and d2 feature vectors
        d1 --- N by 160
        d2 --- N by 160
    return:
        the paired features' index relative to d1 & d2 --- [, , ,] & [ , , ,]
        
    """
    index_source = []
    index_target = []
    tempdis = []
    for pi in range(d1.shape[0]):
        t_p1 = d1[pi,:]
        ds = 10**9
        index_target.append(pi)
        for ti in range(d2.shape[0]):
            temp = sqrt(sum((t_p1-d2[ti,:])**2))
            if temp < ds:
                ds = temp
                target_index = ti
        index_source.append(target_index) # must be at outside of loop 
        tempdis.append(temp)
        
    return index_source, index_target
    
def computePairs(out, out2, index_source, index_target, cDist = 18):
    """
    This function elim
    input:
        paired list of index
        control feature distance
        out ---  feature point outputed by KeyPointsGenerator.py
        index_source --- index of matched pairs
    ouput:
        paired list of points
        
    """
    pairs = []
    vectorDis = []
    
    for ri in range(len(index_target)):
        x, y, z = [], [], []
        x.append(out[0][index_target[ri]])
        x.append(out2[0][index_source[ri]])
    
        y.append(out[1][index_target[ri]])
        y.append(out2[1][index_source[ri]])
        
        z.append(out[2][index_target[ri]])
        z.append(out2[2][index_source[ri]])
        
        x = np.array(x).astype(float)
        y = np.array(y).astype(float)
        z = np.array(z).astype(float)
        
        dist = np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2 + (z[1]-z[0])**2)
        if dist < cDist:
            
            vectorDis.append([x[0]-x[1], y[0]-y[1], z[0]-z[1]])
            
            pairs.append([[x[0],y[0],z[0]],[x[1],y[1],z[1]]])
            
            plot3d(x, y, z, tube_radius=0.25, colormap='Spectral')
            
    arrayPairs = np.array(pairs)
    arrayVector = np.array(vectorDis)
    medianVector = np.median(vectorDis, axis=0)
    
    def removeOutlier(vector, medianVector, crit = 0.9):
        """
        eliminating the outlier direction vectors using rank filter
        """
        mr = np.sqrt(np.sum(medianVector**2))
        medianVector = medianVector/mr
        
        abors = np.sqrt(np.sum(arrayVector**2, axis=1))
        
        tempVector = np.zeros((len(abors),1)) # for the line 87 array division
        tempVector[:,0] = abors # for the line 87 array division
        mrVector = arrayVector/tempVector 
        
        judges = np.dot(mrVector,medianVector)
        
        return judges>crit
    
    boolTruePair = removeOutlier(arrayVector, medianVector)
    
    indexTruePair = np.where(boolTruePair==False)
    
    TruePair = np.delete(arrayPairs,indexTruePair,0)
    
    return TruePair

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
            
def getDisplacements(ref_point, points_arr1, n_points=3):
    """
    ------------------
    this function returns the displacements array calculated 
    using the SVDs of nearest input n_points
    ------------------
    inputs:
        ref_point --- n by 3
        points_arr1 --- n by 3
    outputs:
        dispx --- n by points_ar1.shape[0]
        displacement arrays: [dispx, dispy, dispz]
    
    """
    
    tree = KDTree(ref_point, leaf_size=2)       
    x_disp, y_disp, z_disp= [], [], []
    
    for i in range(point_arr1.shape[0]):
        dist, ind = tree.query([point_arr1[i,:]], k=n_points)
        target = ref_point[ind[0]]
        source = def_point[ind[0]]
        s,r, t = fixed_SVD(target, source)
        x_disp.append(t[0][0])
        y_disp.append(t[0][1])
        z_disp.append(t[0][2])
  
    return [xdisp, ydisp, zdisp]

def showResults(point_arr1, disps, comp="x", geometry=True, image=None):
    """
    showing the computed results
    input: 
        point_arr1 ---  n by 3 array
        disps --- [[],[],[]] list
        comp --- string indicating the component
        geometry --- showing the contour or not
    output:
        None  
    """
    if comp == "x":
        disp = disps[0]
    elif comp == "y":
        disp = disps[1]
    elif comp == "z":
        disp = disps[2]
    else:
        raise TypeError("un-identified input string characters! must be from 'x' , 'y' or 'z'. ")
    
    mesh2 = SpatialSampling.MeshDelaunay(point_arr1.T) # use the original point_arr1 
    mesh2.mesh()
    mesh2.alpha_shape(2*int(Sample_Density))
    mesh2.mesh_triangle() # find unique triangles
        
    for dummy_i in range(1):   
        Mesh_optimizer = SpatialSampling.MeshOptimizer(point_arr1.T,mesh2.ntri)
        Mesh_optimizer.edge_connector_smooth()
        
    if not geometry:
        @mlab.show
        def main():
            mesh2.view(mesh2.polydata(mesh2.points.T,disp))  
        main()
        mayavi.mlab.colorbar(object=None, title="U(Z)(mm)", orientation="vertical", nb_labels=20, nb_colors=None, label_fmt=None)
        "colorbar reference: http://docs.enthought.com/mayavi/mayavi/auto/mlab_decorations.html"    
    
    else:
        if image:
            plot = Visualization.DataVisulization(image, image.mean())
            plot.contour3d()
            surf = plot.surf_points
            inter_z = interpolate.LinearNDInterpolator(mesh2.points, disp, fill_value=0.0, rescale=True)        
            extrapo_z = inter_z.__call__(surf)
            plot.plot3dpoint(extrapo_z)
        else:
            print "need to input image to show the contour!"
            
if __name__ == "__main__":
    print "this is an unit test!"
    
    # generating out,out2
    
    pair_index = minFeaturePairing(d1, d2)
    pairs = computePairs(out, out2, pair_index[0], pair_index[1])

