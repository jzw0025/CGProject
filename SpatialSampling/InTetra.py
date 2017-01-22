"""
This module location each point in the spatial tetrahedrons, and then it computes 
the strain in the element.

"""
import numpy as np

def SameSides(v1,v2,v3,v4, point):
    """
    this function calculates the point if it's in the sameside of each tetrahedorn face
    
    """
    vector1 = np.array(v2) - np.array(v1)
    vector2 = np.array(v3) - np.array(v1)
    normal = np.cross(vector1, vector2)
    
    dot_V4 = np.dot(normal, np.array(v4) - np.array(v1))
    dot_P = np.dot(normal, np.array(point) - np.array(v1))
    #print abs(dot_V4)
    #print dot_V4
    #print dot_V4/abs(dot_V4)
    #print abs(dot_P)
    #print dot_P
    #print dot_P/abs(dot_P)
    return np.sign(dot_V4) == np.sign(dot_P)
    
def inTetra(element_coordinator, point):
    
    """
    the input element_coordinator is a list which each element in the list is a coordinator of a point
    in each dimension 
    
    element_coordinator is 2D list, each element has coordinator in 3D dimensional coordinator system
    
    """
    v1 = element_coordinator[0]
    v2 = element_coordinator[1]
    v3 = element_coordinator[2]
    v4 = element_coordinator[3]
    
    return SameSides(v1,v2,v3,v4,point) & SameSides(v2,v3,v4,v1,point) & SameSides(v3,v4,v1,v2,point) & SameSides(v4,v1,v2,v3,point)
    
    
if __name__ == "__main__":
    """
    building testing cases.
    """
    element_coordinator = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
    point1 = [0.1,0.1,0.1]
    print inTetra(element_coordinator, point1) == True
    point1 = [0.1,0.1,0.0]
    print inTetra(element_coordinator, point1) == False
    point1 = [0.0,0.1,0.1]
    print inTetra(element_coordinator, point1) == False
    point1 = [0.1,0.0,0.1]
    print inTetra(element_coordinator, point1) == False
    point1 = [0.1,0.1,0.3]
    print inTetra(element_coordinator, point1) == True
    point1 = [0.1,0.3,0.1]
    print inTetra(element_coordinator, point1) == True
    point1 = [0.5,0.5,0.5]
    print inTetra(element_coordinator, point1) == False
    point1 = [0.0,1.0,0.0]
    print inTetra(element_coordinator, point1) == False
    point1 = [1.0,0.0,0.0]  
    print inTetra(element_coordinator, point1) == False
    point1 = [0.0,0.0,1.0]
    print inTetra(element_coordinator, point1) == False
    point1 = [100.0,200.0,300.0]
    print inTetra(element_coordinator, point1) == False
    point1 = [-100.0,200.0,300.0]
    print inTetra(element_coordinator, point1) == False
    point1 = [100.0,-200.0,300.0]
    print inTetra(element_coordinator, point1) == False
    point1 = [100.0,200.0,-300.0]
    print inTetra(element_coordinator, point1) == False