'''This code will take any 1D array and convert it into a KDE representation
    Parameters:
        - specific folder containing the 1D lists of choice (folder)
        - covariance factor (covfac)
        - grid spacing for the returned density (delta)
'''

def kernelDensity(Lc, covfac=8, delta=0.25, max_Lc = 600):
    ''' Calculates the kernel density, returns xs and density_xs
    Default covariance factor = 6, delta = 0.5 nm.'''
    from scipy.stats import gaussian_kde
    from numpy import zeros, arange, asarray
    xs = arange(0,max_Lc,delta)
    if (len(Lc) <= 15):
        density_xs = zeros(len(xs))
    else:
        density = gaussian_kde(Lc)
        density.covariance_factor = lambda : covfac/max(Lc)
        density._compute_covariance()
        density_xs = asarray([round(x,10) for x in density(xs)])
    return [xs, density_xs]

def densityList(folder, covfac = 8, delta = 0.25, max_Lc = 600):
    from glob import glob
    list_files = glob("%sList_text/*.txt" % folder)

    density_list = []
    curve_num_list = []
    for file in list_files:
        Lc = load1Col(file)
        density_list.append(kernelDensity(Lc, covfac, delta, max_Lc))
        curve_num_list.append(curveNumFinder(file))

    return density_list, curve_num_list

def saveKDE(folder, covfac = 8, delta = 0.25, max_Lc = 600):
    '''New density folder assumes list folder begins with List_'''
    import os
    density_folder = "%sDensity_text/" % folder
    print density_folder
    if not os.path.isdir(density_folder):
        os.mkdir(density_folder)

    density_list, curve_num_list = densityList(folder, covfac, delta, max_Lc)
    i = 0
    while i < len(density_list):
        xs = density_list[i][0]
        density_xs = density_list[i][1]
        curve = curve_num_list[i]
        filename = "%sDensity_%s.txt" %(density_folder, curve)
        saveFile2(filename, xs, density_xs, 'Value', 'KDE')
        i = i+1

def curveNumFinder(fileName):
    ''' Takes a filename and assuming XXX.txt format, 
        returns the curve number XXX as a string.
        Aborts and returns error if XXX is not a number.'''
    import re
    curve = re.search('[a-z]?[0-9]{3}(?=.txt)', fileName)
    curveNum = curve.group(0)
        
    return curveNum

def load1Col(fileName):
    ''' Takes a file with one columns and returns a list (not an array). '''
    A = []
    file = open(fileName, 'r')
    for line in file:
        if line.startswith('F'):
                continue
        A.append(float(line.split()[0]))
    file.close
    return A

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