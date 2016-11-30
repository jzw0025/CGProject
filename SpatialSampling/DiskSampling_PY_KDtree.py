"""
This module is the pure python implementation for Disksampling module.
It is the copy version for fortran file Sampling.90 file
"""
import random
import numpy as np
from math import *
from scipy import spatial


def distance(tsample, query_point,rck,sample_number,dim,Xlim):
    """
    Inputs:
        
    'tsample' is the whole domain samples
    
    'query_point' (numpy array) is the one sample query point for generating random seeds around it and check the collision
    
    'rck' is the collision distance for checking the samples' space
    
    'sample_number' is the total numbers of random seeds generated around the query point
    
    'dim' is dimensionality of the point
    
    'Xlim' ([[xlow,xhigh],[ylow,yhigh],[zlow,zhigh]]) is the spacial domain constrain for the point generation
    
    Outputs:
        
    'npout' is the new samples would be put into the sample array
    
    'sumiden' is number of new samples (size of none zeros in np)
    
    """

    # for loop: for i in 0,1,2 ... length(sample_number)
    ranp = np.zeros(3)
    check_distance_squared = rck * rck
    identify = np.zeros(sample_number)
    npout = np.zeros([sample_number,dim])
    
    targetTree = spatial.cKDTree(tsample)
    
    rad_points = np.zeros([sample_number,3])
    
    for i in range(sample_number-1,-1,-1):
        # the range third element is the increment!
        # generation random number from [0 1)
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        
        #(1) doing a generation of random number around the fixed query_point
        radius = rck*(r1+1)
        angle1 = 2*pi*r2
        angle2 = 2*pi*r3
        
        # generating x,y,z random coordinates
        #print ranp[0]
        #print query_point[0][0]
        #print radius*cos(angle1)*sin(angle2)
        
        # using 'rck' istead of radius may be more efficient
         
        ranp[0] = query_point[0][0] + rck*cos(angle1)*sin(angle2)
        ranp[1] = query_point[0][1] + rck *sin(angle1)*sin(angle2)
	ranp[2] = query_point[0][2] + rck *cos(angle2)
        
        #ranp[0] = query_point[0][0] + radius*cos(angle1)*sin(angle2)
        #ranp[1] = query_point[0][1] + radius *sin(angle1)*sin(angle2)
	#ranp[2] = query_point[0][2] + radius *cos(angle2)
    
        #(2) check the location of newly generated random number to see whether it's within the domain range
        #  the imported Xlim data formate
        #  [[xlow,xhigh],[ylow,yhigh],[zlow,zhigh]]
        if ranp[0] < Xlim[0,0]:
            ranp[0] = Xlim[0,0]
            
        if ranp[0] > Xlim[0,1]:
            ranp[0] = Xlim[0,1]
            
        if ranp[1] < Xlim[1,0]:
            ranp[1] = Xlim[1,0]
        
        if ranp[1] > Xlim[1,1]:
            ranp[1] = Xlim[1,1]
        
        if ranp[2] < Xlim[2,0]:
            ranp[2] = Xlim[2,0]
        
        if ranp[2] > Xlim[2,1]:
            ranp[2] = Xlim[2,1]
            
        rad_points[i,:] = ranp
            
        #(3) check the collision: checking each point in the samplea array to see if it's collided with the newly generated random point.
        # tsample is transpose, row major 
        
    results = targetTree.query(rad_points) # KdTree Implementation 

    rij =  tsample[results[1]]-rad_points
    
    d2 = np.sum(rij**2, axis=1) # calculate Euclidean Distance  
    identify= np.zeros(sample_number)        
    identify[d2 < check_distance_squared] = 1
    
    # if there is no collision, add it into temporay output array np, and we will add this new point into the existance sample array.    
    # update the number of new sample in the temporary array, adding +1 into the array (suminde = rnum - indentified)
    
    npout =  rad_points[identify==0]     
    # end
             
    sumiden =  sample_number - sum(identify)
    return npout, sumiden

