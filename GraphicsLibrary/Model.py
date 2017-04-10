"""
batch/ vertex list control

"""
import pyglet
import numpy as np
from pyglet.gl import *
import Object

class Model(object):
    
    def __init__(self):
        
        self.batch = pyglet.graphics.Batch()
        self.vertex_list = {} # displaying list
        self._vertex_list = {} # full list
        
    def addObject(self, input_object):
        """
        Object --- input DynamicObjects
        
        """
        if not isinstance(input_object, Object.DynamicObjects):
            raise TypeError("the input variable must be a instance from class Object_R2.DynamicObjects")
        
        def Triangles(glVertices, glIndices, glColor, glNormals=None): 
            """
            help method
            """
            if glNormals is None:
                print "surface normals are not engaged!" 
                #print len(colors)
                #position = input_object.getPosition()
                self._vertex_list[input_object.getPosition()] = self.batch.add_indexed(np.size(glVertices)/3, 
                                                        GL_TRIANGLES,
                                                        None,
                                                        glIndices,
                                                        ('v3f/static', glVertices),                 
                                                        ('c3B/static', glColor))
            else:
                print "found input surface normals!"
                print "indices shape" + str(len(glIndices))
                print "vertices shape" +str(len(glVertices))
                #vertice_normals = list(np.reshape(glNormals.T,np.size(normals.T), order='C'))
                print "normal shape" + str(len(glNormals))
                print "color shape" + str(len(glColor))
                print "location is:" + str(input_object.getPosition())
                self._vertex_list[input_object.getPosition()] = self.batch.add_indexed(np.size(glVertices)/3, 
                                                        GL_TRIANGLES,
                                                        None,
                                                        glIndices,
                                                        ('v3f/static', glVertices),
                                                        ('n3f/static', glNormals),
                                                        ('c3B/static', glColor))
                                                        
        print "found cropped input glVertices list: " + str(len(glVertices))
        print "found cropped input glIndices list: " + str(len(glIndices))
        print "found cropped input glNormals list: " + str(len(glNormals))
        print "found cropped input glColor list: " + str(len(glColor))
                                                        
        glVertices = input_object.getGLvertices()    
        glIndices = input_object.getGLindices()                                             
        glColor = input_object.getGLcolors()
        glNormals = input_object.getGLnormals() 
         
        Triangles(glVertices, glIndices, glColor, glNormals)
        
    def removeObject(self, input_object):
        
        if not isinstance(input_object, Object.DynamicObjects):
            raise TypeError("the input variable must be a instance from class Object_R2.DynamicObjects")
        
        self._vertex_list.pop(input_object.getPosition()).delete()
        
    def showObject(self, input_object):
        
        self.vertex_list[input_object.getPosition()] = self._vertex_list.get(input_object.getPosition())
        
    def hideObject(self, input_object):
        
        self.vertex_list.pop(input_object.getPosition()) 
        
                                                                      
                                                        
        