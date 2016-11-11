"""
GLdata
"""

class GLdata:
    """
    this class loads several 
    
    """
    def __init__(self):
        
        dataList = []
        
        
    def addDataFolder(self, strAddress):
        
        mat_contents1 = sio.loadmat('/Users/junchaowei/Desktop/M6_OD_125_C-scan_full_denoise_Surface_Normal/rim_vertices.mat')
        points = mat_contents1['par1'].T
        mat_contents1 = sio.loadmat('/Users/junchaowei/Desktop/M6_OD_125_C-scan_full_denoise_Surface_Normal/rim_elements.mat')
        elements = mat_contents1['par1']
        mat_contents1 = sio.loadmat('/Users/junchaowei/Desktop/M6_OD_125_C-scan_full_denoise_Surface_Normal/rim_normals.mat')
        normals = mat_contents1['par1'].T
                