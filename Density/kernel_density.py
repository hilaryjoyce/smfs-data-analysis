'''This code will take any 1D array and convert it into a KDE representation
    Parameters:
        - specific folder containing the 1D lists of choice (folder)
        - covariance factor (covfac)
        - grid spacing for the returned density (delta)
        - min and max ranges to place the KDE in
        - the defaults are set up for the WLC protein analysis, NOT DNA analysis, where they should be:
            covfac = 0.3 # covariance factor for KDE of plateau region
            delta = 0.001 # grid size (nN) for calculating the KDE for the plateau region
            min_x = -5.0 # minimum force for density (nN)
            max_x = 5.0 # maximum force for density (could be 1 nN but this allows us to detect wonky curves)
'''

def saveKDE(parameter_folder, covfac = 8, delta = 0.25, min_x = 0, max_x = 600):
    '''Folder must be parameter folder.
       New density folder assumes list folder begins with List_'''
    import os
    density_folder = "%sDensity_text_max%g/" % (parameter_folder, max_x)
    if not os.path.isdir(density_folder):
        os.mkdir(density_folder)

    density_list, curve_num_list = densityList(parameter_folder, covfac, delta, min_x, max_x)
    i = 0
    while i < len(density_list):
        xs = density_list[i][0]
        density_xs = density_list[i][1]
        curve = curve_num_list[i]
        filename = "%sDensity_%s.txt" %(density_folder, curve)
        saveFile2(filename, xs, density_xs, 'Value', 'KDE')
        i = i+1
    return density_folder

def kernelDensity(List, covfac=8, delta=0.25, min_x = 0, max_x = 600):
    ''' Calculates the kernel density, returns xs and density_xs
    Default covariance factor = 6, delta = 0.5 nm.'''
    from scipy.stats import gaussian_kde
    from numpy import zeros, arange, asarray
    xs = arange(min_x,max_x,delta)
    if (len(List) <= 15):
        density_xs = zeros(len(xs))
    else:
        density = gaussian_kde(List)
        if covfac < 3:
            density.covariance_factor = lambda : covfac # /max(List) originally but that doesn't work for peeling plots!
        else: 
            density.covariance_factor = lambda : covfac/max(List)
        density._compute_covariance()
        density_xs = density(xs)
        #density_xs = asarray([round(x,10) for x in density(xs)])
    return [xs, density_xs]

def densityList(folder, covfac = 8, delta = 0.25, min_x = 0, max_x = 600):
    from glob import glob
    list_files = glob("%sList_text/*.txt" % folder)

    density_list = []
    curve_num_list = []
    for file in list_files:
        List = load1Col(file)
        density_list.append(kernelDensity(List, covfac, delta, min_x, max_x))
        curve_num_list.append(curveNumFinder(file))

    return density_list, curve_num_list

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

def load1Col(fileName):
    ''' Takes a file with one columns and returns a list (not an array). '''
    A = []
    file = open(fileName, 'r')
    for line in file:
        if line.startswith('# '):
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