import numpy as np
from skimage import measure
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab 

class VolumeSurface:
    """
    this class return several surfaces for the reconstructed volume
    """
    def __init__(self,volume, *levels):
        self.vol = volume
        self.level = levels # tuple input for different level of isosurfaces (tuple)
        print "found the input level of surfaces" + str(self.level)
        self.surface = [] # store the computed surfaces 
         
    def test(self):
        """
        this is the testing method
        """
        print self.level
        
    def MarchCube(self):
        """
        this method calculates the isosurfaces for the input volume
        """
        def march_cube(volume, level):
            """
            help method for the march cube, using the skimage module
            """
            verts, faces = measure.marching_cubes(volume, level, spacing=(1, 1, 1))
            return verts, faces
        
        dictVolume = {}
        for i in range(len(self.level)):
            templist = []
            print "processing for the isosurface: " + str(self.level[i])
            verts, faces = march_cube(self.vol, self.level[i])
            templist.append(verts)
            templist.append(faces)
            dictVolume[self.level[i]] = templist # output a dictionary with the keys of volume and content of list of surface and vertices
            
        self.surface = dictVolume
        
    def getSurfaces(self):
        return self.surface
        
    def showSurfaces(self,level):
        
        if self.surface == []:
            print "no surfaces are computed in the domain volume!"
            
        print self.surface.keys()    
        
        if level not in self.surface.keys():
            print "could not find the computed surface for the Iso-value: " + str(level)
            return None
        
        surf_lists = self.surface.get(level)
        print surf_lists
        
        # surf_lists[0] --- vertices
        # surf_lists[1] --- surface elements
        verts = surf_lists[0]
        faces = surf_lists[1]
        
        mlab.triangular_mesh([vert[0] for vert in verts],
                             [vert[1] for vert in verts],
                             [vert[2] for vert in verts],
                             faces) 
        mlab.show() 

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],
#                cmap='Spectral', lw=1)
#plt.show()


if __name__ == "__main__":
    print "this is a test!"
    a=[1]
    testclass = VolumeSurface(image1,6000,11000)
    testclass.MarchCube()
    testclass.showSurfaces(11000)