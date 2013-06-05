# Loading essentials

def load2Col(fileName, header=True, col1_num=True):
    from numpy import asarray
    ''' Takes a file with two columns and return each as an array. '''
    A = []
    B = []
    file = open(fileName, 'r')
    if header:
        file.readline()
    i =0
    if col1_num:
        for line in file:
            A.append(float(line.split()[0]))
            B.append(float(line.split()[1]))
        file.close
    else:
        for line in file:
            A.append(line.split()[0])
            B.append(float(line.split()[1]))
        file.close
    A = asarray(A)
    B = asarray(B)
    return [A, B]

def load1Col(fileName, header=True):
    ''' Takes a file with one columns and returns a list (not an array). '''
    A = []
    file = open(fileName, 'r')
    if header:
        file.readline()
    for line in file:
        A.append(float(line.split()[0]))
    file.close
    return A

def loadFolder(folder, columns, header=True):
    from glob import glob
    fileName = folder+"*.txt"
    files = glob(fileName)
    data = []
    if columns == 1:
        for file in files:
            data.append(load1Col(file, header))
    if columns == 2:
        for file in files:
            data.append(load2Col(file, header))
    return data


# Saving essentials

def saveFile1(path, col1, col1_name):
    '''Saves a file with 1 column of data.'''
    file = open(path, 'w')
    file.write("%s\n" % col1_name)
    k  =0
    while k<len(col1):
        el = '%g\n' % col1[k]
        file.write(el)
        k=k+1
    file.close

def saveFile2(path, col1, col2, col1_name, col2_name):
    '''Saves a file with 2 columns of data.'''
    file = open(path, 'w')
    file.write("%s\t%s\n" % (col1_name, col2_name))
    k=0
    while k<len(col1):
        el = '%3g\t%3g\n' % (col1[k], col2[k])
        file.write(el)
        k=k+1
    file.close

# Curve Number identifiers

def curveNumFinder(fileName):
    ''' Takes a filename and assuming {x}XXX.txt or LineXXXXPointXXXXformat, 
        returns the curve number {x}XXX as a string or XXXX.XXXX as the
        curve identifier.'''
    import re
    line_type = re.search('DNA', fileName)
    if line_type:
        curve =  re.findall('[0-9]{4}', fileName)
        curveNum = "%s.%s" % (curve[0], curve[1])
    else:
        curve = re.search('[a-z]?[0-9]{3}(?=.txt)', fileName)
        curveNum = curve.group(0)
    return curveNum

def findPeaks(x, y, tolerance=0.1):
    ''' Takes x and y and find peaks
        that are higher than the tolerance*top peak.
        Returns x, y lists.'''
    import numpy as np
    min_y = max(y)*tolerance
    gradients = np.diff(y)
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
        maxima_x.append(x[index])
        
    maxima_y = []
    for index in max_locations:
        maxima_y.append(y[index])
    
    return maxima_x, maxima_y

def allWLC(Lc):
    '''Calculate worm-like curves for each Lc peak. 
    Returns tss_list, WLC_list, lists of arrays.'''
    WLC_list = []
    tss_list = []
    
    if (len(Lc) == 0): 
        WLC_list.append(0)
    elif (len(Lc) == 1):
        WLC_list.append(calcWLC(Lc))
    else:
        for L in Lc:
            WLC_list.append(calcWLC(L))

    return WLC_list

def calcWLC(Lc, p=0.4):
    ''' Takes the tss array tss and a specific contour length and
        calculates the theoretical force array and array length'''
    from numpy import arange, zeros
    alpha = 1.3806488*10**(-23)*298./(p*(10**(-9))) # (KbT/p, units of m)
    tss = arange(0, 400, 0.1)
    #print Lc
    #print tss[100:110]
    Force = zeros(len(tss))
    i = 0
    while ((tss[i]<Lc) & (i+1<len(tss))):
        Force[i] = alpha*(0.25*(1.-tss[i]/Lc)**(-2) + tss[i]/Lc - 0.25)*10**12
        i+=1
    Force = [x for x in Force[0:i]]
    tss = [x for x in tss[0:i]]
    return [tss, Force]