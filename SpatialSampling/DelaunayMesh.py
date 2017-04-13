""" 
This module contain the delaunay mesh for the 3D sampling points
using dir() to get the module classes included
using .__name__ to find the orignal name of the module if it has been renamed

"""
import numpy as np
from scipy.spatial import Delaunay
import scipy.io as sio
#import DataVisu
#import PoissonDisk  

from numpy import array, random, linspace, pi, ravel, cos, sin, empty
from tvtk.api import tvtk
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi import mlab
       
class MeshDelaunay:
    """
    This class creates tetrahedron element in 3D space, and remove the irregular element using alpha shape function 
    """
    def __init__(self,pointarr):
        self.points = pointarr.T
        # This creates the transpose of point array, the point array dimension is (3,N), which is column major
        # and the transpose of point array dimension is (N,3),which is row major
        self.tri = []
        # this stores the original tetra elements
        self.triangle= []
        # this stores all the triangle facets of tetras
        self.ntri = []
        # after alpha shape operations, the triangles are in ntri
        
    def mesh(self):
        self.tri = Delaunay(self.points)
        # this creates delaunay triangles, and in the 3D space it is a tetrahedron
        # self.tri.simplices.shape is (N,4), which has N elements and 4 vertices
    
    def mesh_triangle(self):
        """
        this method finds the unique triangles from  the delaunay mesh, and then stores those 
        triangles into a self.triangle for data visulization
        """
        print self.ntri.shape
        triangle1 = np.sort(np.array([self.ntri[:,0],self.ntri[:,1],self.ntri[:,2]]), axis = 0)
        triangle2 = np.sort(np.array([self.ntri[:,0],self.ntri[:,1],self.ntri[:,3]]), axis = 0)
        triangle3 = np.sort(np.array([self.ntri[:,0],self.ntri[:,2],self.ntri[:,3]]), axis = 0)
        triangle4 = np.sort(np.array([self.ntri[:,1],self.ntri[:,2],self.ntri[:,3]]), axis = 0)
        # choosing three point to construct triangle from 4 points vertices
        # there are four combinations C(4,3)
        #self.triangle = np.concatenate((triangle1,triangle2,triangle3,triangle4),axis=1)
        # this is tetrahedron element's triangle facets, and it may overlap each other  
        
        triangle_dup = np.concatenate((triangle1,triangle2,triangle3,triangle4),axis=1)
        trans_triangle = triangle_dup.T
        c = np.ascontiguousarray(trans_triangle).view(np.dtype((np.void, trans_triangle.dtype.itemsize * trans_triangle.shape[1])))
        _, idx2 = np.unique(c, return_index=True)
        self.triangle = trans_triangle[idx2]
      
    def alpha_shape(self,radius):
        """
        this method using the alpha shape concepts to find the radius of circumcircle of tetrahedron,
        and it removes the element with R > radius.
        
        """
        simplices = self.tri.simplices
        # Example using index for permute the numpy column
        #>>> i = np.array([0, 2, 4, 3, 1])
        #>>> a 
        #array([[10, 20, 30, 40, 50],
               #[ 6,  7,  8,  9, 10]])
        #>>> a[:, i]
        #array([[10, 30, 50, 40, 20],
               #[ 6,  8, 10,  9,  7]])
        # solve the linear matrix for the center of sphere given four points of tetrahedron
        index = []
        # the following looper is doing the check for each tetra element
        # if the radius of its circumsphere is bigger than the threshold, it will be removed.
        for i in range(self.tri.simplices.shape[0]):
            point1 = self.points[self.tri.simplices[i,0],:]
            point2 = self.points[self.tri.simplices[i,1],:]
            point3 = self.points[self.tri.simplices[i,2],:]
            point4 = self.points[self.tri.simplices[i,3],:]
            AB_vector = point1 - point2
            BC_vector = point2 - point3
            CD_vector = point3 - point4         
            LHS = np.array([AB_vector,BC_vector,CD_vector])
            AB_sqare = sum(point1**2 - point2**2)/2
            BC_sqare = sum(point2**2 - point3**2)/2
            CD_sqare = sum(point3**2 - point4**2)/2
            RHS = np.array([AB_sqare,BC_sqare,CD_sqare])
            center = np.linalg.solve(LHS, RHS)
            # there is referece in the folder for find the Radius of circumsphere of tetrahedron
            R = np.sqrt(sum((center-point1)**2))
            if R > radius:
                index.append(i)
        self.ntri = np.delete(self.tri.simplices,index,axis=0)
            # new tri element
            # print np.allclose(np.dot(LHS, x), RHS)
            # check that the solution is corret.
            
    def polydata(self,point_arr,color=None):
        """
        This method performs a creation of data structure for visulization of tetrahedron.
        It uses the triangle elements for representing the tetrahedron facets
        
        """
        # The numpy array data.
        points = point_arr.T
        triangles = self.triangle # there is transpose, be careful!!!
        #triangles = triangle_set # input variable
        if color is None:
            # default the value for the color, otherwise the color array will be the sample as transpose of input array 'point_arr'
            scalars = random.random(points.shape)
        scalars = color    
        # the scalars is representing the color of the face color of tetrahedrons
        
        # The TVTK dataset.
        mesh = tvtk.PolyData(points=points, polys=triangles)
        mesh.point_data.scalars = scalars
        mesh.point_data.scalars.name = 'scalars'
        return mesh
       
    def view(self,dataset):
        """ Open up a mayavi scene and display the dataset in it.
        """
        fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0),
                      figure=dataset.class_name[3:])
        surf = mlab.pipeline.surface(dataset, opacity=1)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                            color=(0, 0, 0), )
                            
class MeshOptimizer:
    """
    This class deals with the mesh optimization, it applies the concepts of diffusion process(Laplacian Operator)
    The documentation could be found here: 
        (1)http://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/06_smoothing.pdf
        (2)http://w.multires.caltech.edu/pubs/ImplicitFairing.pdf
        (3)http://gustavo.rabello.org/media/pdf/publications/phd/phd2012.pdf
    """
    def __init__(self,pointarr,tri):
        """
        The import point array is in the dimension of (3,N)
        The point array is the transpose of point array
        """
        self.points = pointarr.T
        self.tri = tri 
        # tri is the delaunay triangulation objector
        self.edge = []
        
    def eucli_distance(self,target_point,radius):
        """
        This method tests the distance between the target point and the point cloud.
        Then, it will find out the points within one ring distance.
        """
        import numpy as np
        # points dimension is (N,3)
        number_points = self.points.shape[0]
        # print number_points
        sqare_radius = radius**2
        # the dimension of neighbor point is (3,N)
        neibor_points = np.array([[],[],[]])
        # using sqare distance to check is more efficient.
        for i in range(number_points):
            # target point dimension is (1,3)
            # point cloud dimension is (N,3)
            sqare_dis = sum((target_point - self.points[i,:])**2)
            # find the target point and add it into the list
            if sqare_dis <= sqare_radius and sqare_dis!=0:
                add_point = np.array([[self.points[i,0]],[self.points[i,1]],[self.points[i,2]]])
                neibor_points = np.append(neibor_points,add_point,axis=1)
        return neibor_points
        
    def move_point(self,target_point,neibor_points):
        # the neibor points is dimension (3,N)
        # be careful about the axis of numpy array, they add and subtract results will be different
        # for more reference about np.newaxis and advancing slice 
        # http://docs.scipy.org/doc/numpy/reference/arrays.indexing.html
        xsize,ysize =  neibor_points.shape
        target_point_new = target_point[:, np.newaxis]
        target_point += np.sum(neibor_points - target_point_new, axis=1)
        return target_point
        # move the point according to its neighbors
        
    def mesh_smoother(self):
        number_points = self.points.shape[0]
        radius = 6
        for i in range(number_points):
            print i, number_points
            # print self.points[i,:]
            # the neibors import is (3,N)
            #neibors = Mesh_optimizer.eucli_distance(self.points[i,:],radius)
            neibors = self.eucli_distance(self.points[i,:],radius)
            self.points[i,:] = self.move_point(self.points[i,:],neibors)
            #print self.points[i,:]
    
    def edge_connector_smooth(self):
        
        edge1 = np.sort(np.array([self.tri[:,0],self.tri[:,1]]),axis = 0)
        edge2 = np.sort(np.array([self.tri[:,0],self.tri[:,2]]),axis = 0)
        edge3 = np.sort(np.array([self.tri[:,0],self.tri[:,3]]),axis = 0)
        edge4 = np.sort(np.array([self.tri[:,1],self.tri[:,2]]),axis = 0)
        edge5 = np.sort(np.array([self.tri[:,1],self.tri[:,3]]),axis = 0)
        edge6 = np.sort(np.array([self.tri[:,2],self.tri[:,3]]),axis = 0)
        edge = np.concatenate((edge1,edge2,edge3,edge4,edge5,edge6),axis=1)
        # removing the duplicated edges by first sorting the edges and then remove
        # find the unique reference here: http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        trans_edge = edge.T
        b = np.ascontiguousarray(trans_edge).view(np.dtype((np.void, trans_edge.dtype.itemsize * trans_edge.shape[1])))
        _, idx = np.unique(b, return_index=True)
        self.edge = trans_edge[idx] # this is the unique edges
        
        scalars = random.random(self.points.shape)
        #scalars = np.ones_like(points.shape)*0.5
        # the scalars is representing the color of the face color of tetrahedrons   
        # this is tetrahedron element's triangle facets, and it may overlap each other                                                            
        dt = 1
        lamida = 0.5           
        target_points = self.points  # target_points shold be fixed at each time                                                                                                                                                  
        for i in range(self.points.shape[0]):
            #print i
            connector1 = self.edge[self.edge[:,0] == i,1] # find the connectivity from the first column of edge
            connector2 = self.edge[self.edge[:,1] == i,0] # find the connectivity from the second column of edge
            #if connector1.size == 0 or connector2.size == 0:
            #    break
            #print type(connector1),type(connector2)
            #connector11 = np.asarray(connector1)
            #connector22 = np.asarray(connector2)
            
            #print type(connector11),type(connector22) 
            if connector1.size != 0 and connector2.size !=0:
                connector = np.concatenate((connector1,connector2))
            elif connector1.size == 0:
                connector = connector2
            elif connector2.size == 0:
                connector = connector1
            
            connecting_points = target_points[connector,:]
            
            if connector.shape[0] != 0:
                ##############  uniform laplacian #############
                vector = np.sum(connecting_points - target_points[i,:], axis = 0)/connector.shape[0]
                ############## scale laplacian ################
                #vector = connecting_points - target_points[i,:]
                #norm = np.sqrt(np.sum((connecting_points - target_points[i:])**2, axis = 1))
                # here needs a numpy array broadcasting trick, find a reference here:
                # http://stackoverflow.com/questions/19602187/numpy-divide-each-row-by-a-vector-element
                # it needs np.newaxis or [:,none]
                #norm_vector = vector/norm[:,None]
                #sum_norm_vector = np.sum(norm_vector, axis = 0)
                ############### maybe curvature laplacian ? ##################
                
                ############## chanllegen for the future ############
                self.points[i,:] = target_points[i,:] + lamida * dt * vector  # updating the target points
                # renew the self.points after updating all target_points                
                         
if __name__ == '__main__': 
    """
    This is DelaunayMesh module
    build in module unit test
    """ 
    print "This is a test for tetrahedron mesh." 
    
    import scipy.io as sio
    from scipy import interpolate
    import Visulization
    from scipy import ndimage
    from scipy.interpolate import griddata
    import numpy as np
    import SpatialSampling
    from mayavi import mlab
    import Visulization
    import InputOutput
    
    PathDicom1 = "directory"
    dicom_volume1 = InputOutput.ReadDicom(PathDicom1)
    dicom_volume1.loadfiles()
    image1 = dicom_volume1.DicArray 

    image_trans = SpatialSampling.MeshSize(image1)
    image_trans.grey2binary(4000)
    image_trans.distance_transform()
    
    Samples = SpatialSampling.DiskSamples(image1,image_trans.newimage,4000)
    
    Samples.GenSamples(5,20000)
    point_arr = Samples.samples
    Visulization.DataVisulization(image1,4000).scatterplot(point_arr)
   
    # create discrete samples in the domain using the 4000 greyscale threshold
    mesh2 = MeshDelaunay(point_arr)
    mesh2.mesh()
    mesh2.alpha_shape(5)
    mesh2.mesh_triangle()
    
    print point_arr[:,10]
       
    for dummy_i in range(1):   
       Mesh_optimizer = MeshOptimizer(point_arr,mesh2.ntri)
       Mesh_optimizer.edge_connector_smooth()
        
    print Mesh_optimizer.points.T[:,10]
    # create a 3D tetrahedron mesh for the sampling points
    
    @mlab.show
    def main():
        mesh2.view(mesh2.polydata(Mesh_optimizer.points.T,random.random(point_arr.size/3)))  
    main()
        
                                
      
       
       