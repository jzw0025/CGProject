"""
get clipped region from the clip plane object

"""

import Object
import ClipPlane
import numpy as np

def getVolume(volume, clipPlaneObject, obj_center):
    
    """
    volume  ---  numpy ndarray
    clipPlaneObject --- ClipPlane object
    obj_center --- openGL object center 3 by 1 tuple
    
    """
    
    if not isinstance(clipPlaneObject, ClipPlane.ClipPlane):
        raise TypeError("the input variable must be a instance from class ClipPlane.ClipPlane")
        
    if type(volume) != np.ndarray:
        raise TypeError("the input volume must be numpy array!")
        
    c1 = obj_center[0]
    c2 = obj_center[1]
    c3 = obj_center[2]
    
    x_low = clipPlaneObject.x1 + c1
    x_high = clipPlaneObject.x2 + c1
    y_low = clipPlaneObject.y1 + c2
    y_high = clipPlaneObject.y2 + c2
    z_low = clipPlaneObject.z1 + c3
    z_high = clipPlaneObject.z2 + c3
    
    xl, yl, zl = volume.shape
    
    # make all values valid
    if x_low < 0: x_low =0 
    if x_high > xl: x_high = xl
    if y_low < 0: y_low =0 
    if y_high > xl: y_high = xl
    if z_low < 0: z_low =0 
    if z_high > xl: z_high = xl
    
    volume_index = np.zeros_like(volume)
    
    print x_low, x_high, y_low, y_high, z_low, z_high
    volume_index[x_low:x_high, y_low:y_high, z_low:z_high] = True
    
    print volume.mean()
    print sum(volume_index)
    
    outvolume = volume * volume_index
     
    print (outvolume).mean()
    
    return outvolume
          

if __name__ == "__main__":
    print "this is unit test!"
    import scipy.io as sio
    
    volume1 = deformed1
    volume2 = deformed2
    
    clip_plane = ClipPlane.ClipPlane()
    clip_plane.applyCenterPoint(0,-50,-80)
    obj_center = (326.55084228515625, 301.58139038085938, 199.94027709960938)
    clip_volume1 = getVolume(volume1,clip_plane,obj_center)
    clip_volume2 = getVolume(volume2,clip_plane,obj_center)
    
    