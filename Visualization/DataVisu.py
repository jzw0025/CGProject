import numpy as np
import mayavi.mlab as p3d 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class DataVisulization():
    def __init__(self, volume,iso):
        self.scalars = volume
        self.surf_points = []
        self.surf_elements = []
        self.color = []
        self.iso = iso
        # self.colordatalen = numpy.linspace(-pi,pi,len(poly_data_object.point_data.scalars))

    def contour3d(self):
        #create one isosurface data set 
        # mayavi.mlab.figure(figure=None, bgcolor=None, fgcolor=None, engine=None, size=(400, 350))
        my_obj = p3d.contour3d(self.scalars, contours=[self.iso]) 
        my_actor=my_obj.actor.actors[0] 
        poly_data_object=my_actor.mapper.input 
        #array with all vertices 
        the_points = np.array(poly_data_object.points)
        self.surf_points = the_points
        # (number_of_points by 3) 
        #  poly_data_object.polys.data.to_array() gives a 1D array for cell data 
        #  of the form [n_vertices_in_face0, 
        #             vert0,vert1,vert2,..., 
        #             n_vertices_in_face1,vert0,vert1,...] 
        # If all faces are triangles, then the following works 
        # The first column gives number of vertices for the poly 
        
        # newshape : int or tuple of ints
        # The new shape should be compatible with the original shape. 
        # If an integer, then the result will be a 1-D array of that length. 
        # One shape dimension can be -1. In this case, the value is inferred from 
        # the length of the array and remaining dimensions. e.g. [-1,4]
        
        the_cells = np.reshape(poly_data_object.polys.data.to_array(), [-1,4]) 
        #Since all faces are triangles, get rid of first column 
        the_cells = the_cells[:,1:]
        self.surf_elements = the_cells
    
    def colorinterpolation():
        pass
        
    def plot3dpoint(self,color):
        #print color.shape
        simple_tm=p3d.triangular_mesh(self.surf_points[:,0], self.surf_points[:,1], 
                               self.surf_points[:,2], 
                               self.surf_elements, 
                               scalars=color)                       
        p3d.colorbar(object=simple_tm, title='Displacement', orientation=None, nb_labels=5, nb_colors=25, label_fmt='%.2f')
                               
    def scatterplot(self,point_arr):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xs = point_arr[0]
        ys = point_arr[1]
        zs = point_arr[2]
        ax.scatter(xs, ys, zs, c='blue', marker='o')   
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_aspect('auto')
        plt.show()