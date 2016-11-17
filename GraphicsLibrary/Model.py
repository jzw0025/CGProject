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
        self._vertex_list = {} # full list
        self._vertex_list_track = {} # keep track of properties of vertex list if is full object or cropped
        
    def addObject(self, input_object):
        """
        Object --- input DynamicObjects
        
        """
        
        # check if the object is already in the vertex list, otherwise add the new object(location)
        print input_object.getPosition() 
        print self._vertex_list.keys()
        if input_object.getPosition() in self._vertex_list.keys():
            print "the object is already exist in the vertex list"
            return
        else:    
            print "input a new object!"            
            if not input_object.getState(): # if not the global object
                self.removeGlobalObject() # remove all the global object when there is local object imported

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
                                                        
                self._vertex_list_track[input_object.getPosition()] = input_object.getState()                                  
                                                
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
                                                        
                self._vertex_list_track[input_object.getPosition()] = input_object.getState()                                        
                                                        
        glVertices = input_object.getGLvertices()    
        glIndices = input_object.getGLindices()                                             
        glColor = input_object.getGLcolors()
        glNormal = input_object.getGLnormals() 
        
        print "found cropped input glVertices list: " + str(len(glVertices))
        print "found cropped input glIndices list: " + str(len(glIndices))
        print "found cropped input glNormals list: " + str(len(glNormal))
        print "found cropped input glColor list: " + str(len(glColor))
         
        Triangles(glVertices, glIndices, glColor, glNormals=glNormal)
        
        
    def removeGlobalObject(self):
        print "removing global objects"
        """
        remove all the gloabl object
        
        """
        print self._vertex_list.keys()
        print self._vertex_list_track.keys()
        for i in range(len(self._vertex_list_track.keys())):
            
            print self._vertex_list_track.get(self._vertex_list_track.keys()[i])
            
            if self._vertex_list_track.get(self._vertex_list_track.keys()[i]): # if it is a global object   
                print i, self._vertex_list_track.keys()[i]
                self._vertex_list.pop(self._vertex_list_track.keys()[i]).delete()
                self._vertex_list_track.pop(self._vertex_list_track.keys()[i])    

    def removeObject(self, input_object):
        
        if not isinstance(input_object, Object.DynamicObjects):
            raise TypeError("the input variable must be a instance from class Object_R2.DynamicObjects")
        
        self._vertex_list.pop(input_object.getPosition()).delete()
        
    def showObject(self, input_object):
        
        self.vertex_list[input_object.getPosition()] = self._vertex_list.get(input_object.getPosition())
        
    def hideObject(self, input_object):
        
        self.vertex_list.pop(input_object.getPosition()) 
        
                                                                      
                                                        
        