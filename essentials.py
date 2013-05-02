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

def curveNumFinder(filename):
    ''' Takes a filename and assuming lineXXXXpontXXXX.txt format, 
        returns the curve number XXXX.XXXX as a string.'''
    import re
    curve =  re.findall('[0-9]{4}', filename)
    curveNum = "%s.%s" % (curve[0], curve[1])
    return curveNum