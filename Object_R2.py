"""
this module contains class used in the pyglet dynamic object rendering
"""
import numpy as np
from pyglet.gl import *
import ClipPlane

class Object():
    """
    Object Class
    This class provides the basic object for displaying in the pyglet
    """
    def __init__(self,center, color):
        self.center = center
        self.color = color
        
    def getObjectCenter(self):
        return self.center
        
    def getObjectColor(self):
        return self.color
        
        
class DynamicObjects(object):
    """
    keep track of dynamically drawing objects
    this should include the batch argument 
    quad, surface(triangles)
    
    """
    
    def __init__(self, points, elements, normals=None, color=None, center=None):
        """
        points --- (3,N) numpy array
        element --- (3,N) or (4,N) numpy array
        normals --- (3,N) numpy array
        color --- (3,-) tuple or list
        center --- (3,-) tuple or list 
        
        """
        ## center
        if center is None:
            self.center = [points[0,:].mean(), points[1,:].mean(), points[2,:].mean()]
        else:
            self.center = center
        
        ## color
        if color is None:
            self.color = (255,0,0) # default red color
        else:
            self.color = color
                    
        ###  points ###
        if type(points) != np.ndarray:
            raise TypeError("the input points type is not numpy nd array")
        else:
            if points.shape[0] != 3:
                self.points =  points.T # transpose the points array into (N,3) array
                self.points[0,:] = self.points[0,:] - self.points[0,:].mean() + self.center[0] # centralize the points array
                self.points[1,:] = self.points[1,:] - self.points[1,:].mean() + self.center[1] # apply the new center
                self.points[2,:] = self.points[2,:] - self.points[2,:].mean() + self.center[2]
            else:
                self.points = points
                self.points[0,:] = self.points[0,:] - self.points[0,:].mean() + self.center[0] # centralize the points array
                self.points[1,:] = self.points[1,:] - self.points[1,:].mean() + self.center[1] # apply the new center
                self.points[2,:] = self.points[2,:] - self.points[2,:].mean() + self.center[2]
     
        ### elements ###
        if type(elements) != np.ndarray:
            raise TypeError("the input elements type is not numpy nd array")
        else:
            if elements.shape[0] == 3:
                print "the input element shape is not (N,3), it is: " + str(elements.shape)
                self.elements =  elements.T # transpose the points array into (N,3) array
            elif elements.shape[0] == 4:
                print "the input element shape is not (N,4), it is: " + str(elements.shape)
                self.elements =  self.extractTriangle(elements.T) # transpose the element and extract surface element into (N,3) array 
            else:
                print "import element size: " + str(elements.shape)
                self.elements = elements
                                 
        ### input normals ###  
        self.normals = normals  
        if self.normals is not None:
            if type(self.normals) != np.ndarray:
                raise TypeError("the input normals type is not numpy nd array")  
            elif self.normals.shape != self.points.shape:
                    print "normals' shape:" + str(self.normals.shape)
                    print "points' shape" + str(self.points.shape)
                    raise ValueError("the input normals do not match with the points!")
            else:
                self.gl_normals = list(np.reshape(self.normals.T,np.size(self.normals.T), order='C'))
                              
        self.gl_vertices = list(np.reshape(self.points.T,np.size(self.points.T), order='C')) #indices = [0,1,2,0,2,3,0,3,1,1,2,3]
        self.gl_indices = list(np.reshape(self.elements, np.size(self.elements), order = 'C'))
        color_array = np.ones((self.points.T.shape[0],3))
        color_array[:,0] = self.color[0]
        color_array[:,1] = self.color[1]
        color_array[:,2] = self.color[2]
        int_color_array = color_array.astype(int) 
        self.gl_colors = list(np.reshape(int_color_array, np.size(int_color_array), order = 'C')) 
          
    def getGLvertices(self):
        return self.gl_vertices
        
    def getGLindices(self):
        return self.gl_indices
    
    def getGLcolors(self):
        return self.gl_colors
        
    def getGLnormals(self):
        return self.gl_normals
        
    @staticmethod                
    def extractTriangle(tetras):
        """
        this help method creates triangle elements from the tetrahedron element
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
                                                    
    def applyClipPlane(self, clipPlaneObject):
        
        if not isinstance(clipPlaneObject, ClipPlane.ClipPlane):
            raise TypeError("the input variable must be a instance from class ClipPlane.ClipPlane")
            
        x_low = clipPlaneObject.x1 
        x_high = clipPlaneObject.x2
        y_low = clipPlaneObject.y1
        y_high = clipPlaneObject.y2
        z_low = clipPlaneObject.z1
        z_high = clipPlaneObject.z2
          
        x_axis = np.logical_and(self.points[0,:] >= x_low, self.points[0,:] <= x_high)
        y_axis = np.logical_and(self.points[1,:] >= y_low, self.points[1,:] <= y_high)
        z_axis = np.logical_and(self.points[2,:] >= z_low, self.points[2,:] <= z_high)
        
        #print sum(x_axis), sum(y_axis), sum(z_axis)
        
        clip_logic = np.logical_and(x_axis, y_axis, z_axis)   
        
        index = np.array(range(self.points.shape[1])) # create index for referencing
        
        new_points = self.points[:,clip_logic] # update points
        
        if self.normals is not None:
            new_normals =  self.normals[:,clip_logic] # update normals
            print "the cropped normals: " + str(new_normals.shape)
            print "the cropped points: " + str(new_points.shape)
        
        active_point_number = set(index[clip_logic]) # convert the index array into a set
        
        print "found activated points number: " + str(len(active_point_number))
        new_elements = []
        #print self.elements.shape  
        for i in range(self.elements.shape[0]):
            setElement = set(self.elements[i,:])
            if len(setElement.intersection(active_point_number))>=3:
                new_elements.append(self.elements[i,:])
       
        #    judge = True # False -- do not add this point
        #    for j in range(self.elements.shape[1]):
        #        if self.elements[i,j] not in active_point_number:
        #            judge = False
        #            break              
        #    if judge:
        #        new_elements.append(self.elements[i,:])
        
        new_elements = np.array(new_elements)
        
        print "applied clipped plane!"
        print "new points" + str(new_points.shape) + " old points" + str(self.points.shape)
        print "new elements" + str(new_elements.shape)+ " old elements" + str(self.elements.shape)
        print "new normals" + str(new_normals.shape)+ " old normals" + str(self.normals.shape)
                
        self.points = new_points # update the normals
        self.elements = new_elements
        self.normals = new_normals
        
        ## update the gl components
        self.gl_vertices = list(np.reshape(self.points.T,np.size(self.points.T), order='C')) #indices = [0,1,2,0,2,3,0,3,1,1,2,3]
        self.gl_indices = list(np.reshape(self.elements, np.size(self.elements), order = 'C'))
        color_array = np.ones((self.points.T.shape[0],3))
        color_array[:,0] = self.color[0]
        color_array[:,1] = self.color[1]
        color_array[:,2] = self.color[2]
        int_color_array = color_array.astype(int) 
        self.gl_colors = list(np.reshape(int_color_array, np.size(int_color_array), order = 'C')) 
        self.gl_normals = list(np.reshape(self.normals.T,np.size(self.normals.T), order='C'))
            
if __name__ == "__main__":
    print "this is unit test!" 
    clip_plane = ClipPlane.ClipPlane()
    clip_plane.applyCenterPoint(points[:,0].mean(),points[:,1].mean(),points[:,2].mean()) 
    dobj = DynamicObjects(points, elements,normals=points.T)
    dobj.applyClipPlane(clip_plane)
    dobj.getGLvertices()
    
           
        