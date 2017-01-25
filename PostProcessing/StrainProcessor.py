"""
this module runs the strain and post processing for the 
Digital Volume Correlation
"""
import numpy as np
from scipy import interpolate

def interpDisplacement(points, displacements, sample_length, coordinators):
    """
    this function does the interpolation over the neighbor connecting points
    input:
        points --- (N,3)
        displacements --- (3, N)  List--Array
        sample_length --- float
        coordinators --- (1,3)
    output:
        strains --- (1,N)
    """
    U = interpolate.LinearNDInterpolator(points, displacements, fill_value=0.0)
    V = interpolate.LinearNDInterpolator(points, displacements, fill_value=0.0)
    W = interpolate.LinearNDInterpolator(points, displacements, fill_value=0.0)
    xcor = coordinators[0]
    ycor = coordinators[1]
    zcor = coordinators[2]
    
    return inter_z
    
def calculateStrain(points, ntri, components, num_smoothes):
        """
        input:
            ntri --- input element (N,3)
            num_smoothes --- int N
            components "containing the components needed to be smoothed" ---- List [[np.array],[np.array],...] 
            points -- (N,3)
            
        """
        edge1 = np.sort(np.array([ntri[:,0], ntri[:,1]]),axis = 0)
        edge2 = np.sort(np.array([ntri[:,0], ntri[:,2]]),axis = 0)
        edge3 = np.sort(np.array([ntri[:,0], ntri[:,3]]),axis = 0)
        edge4 = np.sort(np.array([ntri[:,1], ntri[:,2]]),axis = 0)
        edge5 = np.sort(np.array([ntri[:,1], ntri[:,3]]),axis = 0)
        edge6 = np.sort(np.array([ntri[:,2], ntri[:,3]]),axis = 0)
        edge = np.concatenate((edge1,edge2,edge3,edge4,edge5,edge6),axis=1)
        trans_edge = edge.T
        b = np.ascontiguousarray(trans_edge).view(np.dtype((np.void, trans_edge.dtype.itemsize * trans_edge.shape[1])))
        _, idx = np.unique(b, return_index=True)
        mesh_edge = trans_edge[idx] # this is the unique edges
        for nu in range(num_smoothes):
            for i in range(points.shape[0]):
                print i
                connector1 = mesh_edge[mesh_edge[:,0] == i,1] # find the connectivity from the first column of edge
                connector2 = mesh_edge[mesh_edge[:,1] == i,0] # find the connectivity from the second column of edge
                if connector1.size !=0  and connector2.size != 0:
                    connector = np.concatenate((connector1,connector2))
                elif connector1.size == 0:
                    connector = connector2
                elif connector2.size == 0:
                    connector = connector1
                else:
                    print "the point is isolated without any connecting edge!"
                connecting_points = points[connector,:]
                connecting_displacements = displacements[connector,:]
                interpDisplacement(connecting_points, connecting_displacements, sample_length, coordinators)
                
                    
                ## smoothing for each strain elements.
                connecting_s11 = s11[connector]
                connecting_s12 = s12[connector]
                connecting_s13 = s13[connector]
                connecting_s22 = s22[connector]
                connecting_s23 = s23[connector]
                connecting_s33 = s33[connector]
                    
                #print connecting_xdis.shape
                #connecting_ydis = extrapo_y[connector,:]
                #connecting_zdis = extrapo_z[connector,:]
                if connector.shape[0] != 0:
                ##############  uniform smoothing #############
                    s11[i] = (np.sum(connecting_s11)+s11[i])/(connector.shape[0]+1) # plus 'one' is the central displace point 
                    s12[i] = (np.sum(connecting_s12)+s12[i])/(connector.shape[0]+1) 
                    s13[i] = (np.sum(connecting_s13)+s13[i])/(connector.shape[0]+1) 
                    s22[i] = (np.sum(connecting_s22)+s22[i])/(connector.shape[0]+1) 
                    s23[i] = (np.sum(connecting_s23)+s23[i])/(connector.shape[0]+1) 
                    s33[i] = (np.sum(connecting_s33)+s33[i])/(connector.shape[0]+1) 
                    #disp[i] = np.median(connecting_xdis)
                    #if abs(s[i]) >1:
                    #    break
                    

def computeStrains():
    pass

def smoothStrains(strain):
    pass
    