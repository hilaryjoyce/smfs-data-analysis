'''This code will find peaks and convert those into useful forms of information to 
    - locate WLCs 
    - number of peaks
    - force jump heights
'''

def findPeaks(x, y, tolerance=0.08):
    from numpy import diff
    ''' Takes x and y and find peaks that are higher than the .
        Returns x, y lists. This is a lame peak-finding algorithm'''
    min_y = tolerance
    gradients = diff(y)
    maxima_num = 0
    max_locations = []
    count=0
    
    for i in gradients[:-1]:
        count+=1
        if ((cmp(i,0)>0) & (cmp(gradients[count],0)<0) & (i != gradients[count])):
            if (y[count] > min_y):
                maxima_num+=1
                max_locations.append(count)

    maxima_x = []
    for index in max_locations:
        maxima_x.append(round(x[index],9))
        
    maxima_y = []
    for index in max_locations:
        maxima_y.append(round(y[index],9))
    
    return maxima_x, maxima_y

def findPeaks_second(x, y, tolerance = 0.00005):
    '''Finds the second derivative and then the minima to locate all peaks.
        Requires smooth data.'''
    from numpy import gradient
    max_x, max_dd_y = findPeaks(x,-1*gradient(gradient(y)), tolerance)
    return max_x, max_dd_y

def jumps(x, to_zero = True):
    '''Takes an array x and returns the difference between values,
    and between the last value in the array and 0.'''
    jumps = []
    if len(x) == 1:
        jumps.append(0-x[0])
    else:
        i = 1
        while i < len(x):
            jumps.append(x[i] - x[i-1])
            i = i+1
        if to_zero:
            jumps.append(0-x[i-1])
    return jumps