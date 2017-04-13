"""
this module loads several data address, and creates a database

"""

import scipy.io as sio
import numpy as np

class DVCdata():
    """
    this is standard input and output data for the DVC
    The data must be saved according the routine in saveDVC module of SaveLoad
    
    """
    def __init__(self, addressFolder):
        
        try:
            self.pre_image1 = sio.loadmat(addressFolder + "/image1")
            #print "imported dictionary: " + str(self.pre_image1)
            for i in range(len(self.pre_image1 )):
                if type(self.pre_image1.get(self.pre_image1.keys()[i])) == np.ndarray:
                    self.image1 = self.pre_image1.get(self.pre_image1.keys()[i])
        except IOError:
            print "can not find the file named: " + addressFolder + "image1"
        
        try: 
            self.pre_image2 = sio.loadmat(addressFolder + "/image2")
            #print "imported dictionary: " + str(self.pre_image1)
            for i in range(len(self.pre_image2)):
                if type(self.pre_image2.get(self.pre_image2.keys()[i])) == np.ndarray:
                    self.image2 = self.pre_image2.get(self.pre_image2.keys()[i])
                    
        except IOError:
            print "can not find the file named: " + addressFolder + "image2"
        
        try: 
            self.pre_point_arr1 = sio.loadmat(addressFolder + "/point_arr1")
            self.point_arr1 = self.pre_point_arr1[self.pre_point_arr1.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "point_arr1"
        
        try: 
            self.pre_point_arr2 = sio.loadmat(addressFolder + "/point_arr2")
            self.point_arr2 = self.pre_point_arr2[self.pre_point_arr2.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "point_arr2"
            
        try: 
            self.pre_surf1_elements = sio.loadmat(addressFolder + "/surf1_elements")
            self.surf1_elements = self.pre_surf1_elements[self.pre_surf1_elements.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "surf1_elements"
        
        try: 
            self.pre_surf1_points = sio.loadmat(addressFolder + "/surf1_points")
            self.surf1_points = self.pre_surf1_points[self.pre_surf1_points.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "surf1_points"
        
        try: 
            self.pre_surf2_elements = sio.loadmat(addressFolder + "/surf2_elements")
            self.surf2_elements = self.pre_surf2_elements[self.pre_surf2_elements.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "surf2_elements"
        
        try: 
            self.pre_surf2_points = sio.loadmat(addressFolder + "/surf2_points")
            self.surf2_points = self.pre_surf2_points[self.pre_surf2_points.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "surf2_points"
            
        try: 
            self.pre_meshPoints = sio.loadmat(addressFolder + "/mesh_points")
            self.meshPoints = self.pre_meshPoints[self.pre_meshPoints.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "mesh_points"
            
        try: 
            self.pre_meshElements = sio.loadmat(addressFolder + "/mesh_elements")
            self.meshElements = self.pre_meshElements[self.pre_meshElements.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "mesh_elements"
       #
       # the mesh is created according to the point_arr1, which has one step smooth
        try: 
            self.pre_dispx = sio.loadmat(addressFolder + "/dispx")
            self.dispx = self.pre_dispx[self.pre_dispx.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "dispx"
            
        try: 
            self.pre_dispy = sio.loadmat(addressFolder + "/dispy")
            self.dispy = self.pre_dispy[self.pre_dispy.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "dispy"
            
        try: 
            self.pre_dispz = sio.loadmat(addressFolder + "/dispz")
            self.dispz = self.pre_dispz[self.pre_dispz.keys()[2]]
            
        except IOError:
            print "can not find the file named: " + addressFolder + "dispz"
   
    def getImage1(self):
        return self.image1
        
    def getImage2(self):
        return self.image2
        
    def getPoints1(self):
        return self.point_arr1
    
    def getPoints1centered(self):
        point_arr1 = np.zeros_like(self.point_arr1)
        point_arr1[:,0] = self.point_arr1[:,0] - self.point_arr1[:,0].mean()
        point_arr1[:,1] = self.point_arr1[:,1] - self.point_arr1[:,1].mean()
        point_arr1[:,2] = self.point_arr1[:,2] - self.point_arr1[:,2].mean()
        return point_arr1
    
    def getPoints2(self):
        return self.point_arr2
        
    def getSurfPoint1(self):
        return self.surf1_points
    
    def getSurfPoint2(self):
        return self.surf2_points
        
    def getSurfElement1(self):
        return self.surf1_elements
    
    def getSurfElement2(self):
        return self.surf2_elements
        
    def getMeshPoints(self):
        return self.meshPoints
    
    def getMeshElements(self):
        """
        the Elements is tetrahedron
        """
        return self.meshElements 
        
    def getDisplacementX(self):
        return self.dispx
    
    def getDisplacementY(self):
        return self.dispy
    
    def getDisplacementZ(self):
        return self.dispz
           
if __name__ == "__main__":
    print "this is the test!"
    
    testClass = DVCdata("Dir")
        
    print testClass.getImage1()
    
    print type(testClass.getImage2())
    
    print testClass.getPoints1()
    
    print testClass.getPoints2()
    
    print testClass.getMeshElements()
    print "elements"
    
    
        