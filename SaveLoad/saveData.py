"""
this module saves the computed Python-DVC data into a readable DVC database.
The class "DVCdata" loads this data in the folder.

"""
import sys
import os
import scipy.io as sio

class saveDataBase():
    
    def __init__(self, address):
        
        if type(address) != str:
            print "The input type must be a string!"     
        
        else:
            assert os.path.exists(address), "I did not find the directory:, "+str(address)     
            self.directoryFileName = address
             
    def save(self, dict_data):
        """
        the input is given by a dict: {"image1": image1, "image2": image2, "point_arr1": point_arr1, "point_arr2" : point_arr2,
                                        "referennce_point_arr1", referennce_point_arr1, "referennce_point_arr2": referennce_point_arr2,
                                        "dispx":dispx, "dispxx": dispxx,
                                        "dispy":dispy, "dispyy": dispyy,
                                        "dispz":dispz, "dispzz": dispzz,
                                        "surf1_points",surf1_points, "surf2_points",surf2_points,
                                        "surf1_elements",surf1_elements,"surf2_elements",surf2_elements,
                                        "mesh_points",mesh_points, "mesh_elements", mesh_elements}
        
        this method reads all the data from the DVC computation, and save those data into a file.
        Notice: the name space is fixed according the programming.
        
        """
        for i in range(len(dict_data.keys())):
            # save the database according to the dictionaries' name
            print dict_data.get(dict_data.keys()[i])
            sio.savemat(self.directoryFileName + dict_data.keys()[i], {'par1':dict_data.get(dict_data.keys()[i])})
        
        print "saved the data into: " + self.directoryFileName
            
            
if __name__ == "__main__":
    import numpy as np
    test = saveDataBase("/Users/junchaowei/Desktop/test/")
    image1 = np.array([[1],[2]])
    image2 = 2
    point_arr1 = 3
    point_arr2 = 4
    referennce_point_arr1 = 5
    referennce_point_arr2 = 6
    dispx = 7
    dispxx = 8
    dispy = 9
    dispyy = 10
    dispz = 11
    dispzz = 12
    surf1_points = 13
    
    datadict={}
    
    datadict['image1'] = image1
    datadict['image2'] = image2
    datadict['point_arr1'] = point_arr1
    datadict['point_arr2'] = point_arr2
    datadict['referennce_point_arr1'] = referennce_point_arr1
    datadict['referennce_point_arr2'] = referennce_point_arr2
    datadict['dispx'] = dispx
    datadict['dispxx'] = dispxx
    datadict['dispy'] = dispy
    datadict['dispyy'] = dispyy
    datadict['dispz'] = dispz
    datadict['dispzz'] = dispzz
    datadict['surf1_points'] = verts # from mlab march cube
    datadict['surf2_points'] = verts2 # from mlab march cube
    datadict['surf1_elements'] = faces
    datadict['surf2_elements'] = faces2                    
    
    save.save(datadict)
        
        
        





        