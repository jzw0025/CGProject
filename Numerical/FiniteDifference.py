"""
this file contains numerical functions

"""
def calculateFirstDerivative(xdata, ydata):
    """
    this function calculates the derivatives of given input data
    input
        xdata --- [1,N] list
        ydata --- [1,N] list
    output
        ydata_dxx --- [1,N] list
    """
    
    if len(xdata) != len(ydata):
        raise TypeError("the input data dimension mismatch!")
    ydata_dxx = []
    for i  in range(len(ydata)):
        if i == 0:
            dxx = (ydata[i+1] - ydata[i])/(xdata[i+1] - xdata[i])
            ydata_dxx.append(dxx)
        elif i == len(ydata)-1:
            dxx = (ydata[i] - ydata[i-1])/(xdata[i] - xdata[i-1])
            ydata_dxx.append(dxx)
        else:
            dxx = (ydata[i+1] - ydata[i-1])/(xdata[i+1] - xdata[i-1])
            ydata_dxx.append(dxx)
    return ydata_dxx
    
def calculateSecondDerivative(xdata, ydata):  
        """
        this function calculates the derivatives of given input data
        input
            xdata --- [1,N] list
            ydata --- [1,N] list
        output
            ydata_dxx --- [1,N] list
        """
        if len(xdata) != len(ydata):
            raise TypeError("the input data dimension mismatch!")
        ydata_dxx = []
        for i  in range(len(ydata)):
            if i == 0:
                dxx = (ydata[i+1] + ydata[-1] - 2*ydata[i])/(xdata[i+1] - xdata[i])/2
                ydata_dxx.append(dxx)
            elif i == len(ydata)-1:
                dxx = (ydata[i-1] + ydata[0] - 2*ydata[i])/(xdata[i] - xdata[i-1])/2
                ydata_dxx.append(dxx)
            else:
               dxx = (ydata[i-1] + ydata[i+1] - 2*ydata[i])/(xdata[i+1]- xdata[i-1])
               ydata_dxx.append(dxx)
        return ydata_dxx
