import numpy as np

"""
This module calculates the mean value and smoothes the data.

"""
### Running mean/Moving average
def running_mean(l, N):
    sum = 0
    result = list( 0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result

### Neighbor Mean/Moving Average
def neighbor_mean(input_list, N):
    """
    inputs:
        input_list --- [] 
        N --- number of neighbors to be averaged
    outputs:
        result --- []
    """
    if  type(input_list) != list:
        raise TypeError("The input list must be a list!")
        
    if type(N) != int:
        raise TypeError("The input type must be an int!")
        
    result = []
    radius = (N-1)/2

    for i in range(len(input_list)):
  
        if (i-radius)<0:
            result.append(np.mean(input_list[:i+1])) 
            
        elif (i+radius)>len(input_list)-1:
            result.append(np.mean(input_list[i:]))
               
        else:
            temp = []
            for j in range(i-radius, i+radius+1):
                temp.append(input_list[j])
            result.append(np.mean(temp))

    return result
    
if __name__== "__main__":
    print "This is a test!"
    x = [i for i in range(-10,11)]
    a = [i**2 for i in range(-10,11)]
    print a 
    sa = neighbor_mean(a,5)
    print sa
    import matplotlib.pyplot as plt
    plt.figure(0)
    plt.plot(x,a,'r',x,sa,'b')
    plt.show()
    
    
    