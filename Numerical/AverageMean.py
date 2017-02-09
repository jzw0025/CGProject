import numpy as np


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
    print type(input_list)
    result = []
    radius = (N-1)/2

    for i in range(len(input_list)):
        if i == 0:
            result.append((input_list[i]+input_list[i+1])/2) 
        elif i == len(input_list)-1:
            result.append((input_list[i-1]+input_list[i])/2)   
        else:
            temp = []
            for j in range(i-radius, i+radius+1):
                temp.append(input_list[j])
            result.append(np.mean(temp))

    return result