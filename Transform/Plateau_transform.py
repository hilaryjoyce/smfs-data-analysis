'''This code WILL take force plateau data in the force-tss form and 
    transform it into a list of valid forces.
    Old parameters I was using that should come in handy:
        baseline_SD = 0.013 # standard deviation for finding baseline
        plat_covfac = 0.3 # covariance factor for KDE of plateau region
        plat_delta_force = 0.001 # grid size (nN) for calculating the KDE for the plateau region
        plat_peak_minimum = 5.0 # minimum peak height for plateau KDE peak finding
        minimum_force = -5.0 # minimum force for density (nN)
        maximum_force = 5.0 # maximum force for density (could be 1 nN but this allows us to detect wonky curves)

'''

def PlateauForceList(folder, SD = 0.017):
    '''Find the list of force values for the long plateau at the end of the curve.'''
    from glob import glob

    tss_force_files = glob('%sTss_Force_text/*.txt' % folder)

    force_list = []
    curve_num_list = []

    for filename in tss_force_files:
        tss, force = load2Col(filename)
        curveNum = curveNumFinder(filename)
        start, stop = start_stop(tss, force, SD)

        mid = start
        while force[mid] < -1.0:
            mid = mid + 1

        force_list.append(force[mid:stop])
        curve_num_list.append(curveNum)

    return force_list, curve_num_list

def saveForceList(folder, SD = 0.017):
    import os
    from numpy import zeros

    if not os.path.isdir("%sTransform_analysis/" % folder):
        os.mkdir("%sTransform_analysis/" % folder)

    parameter_folder = "Parameters_SD%g" % SD
    if not os.path.isdir("%sTransform_analysis/%s" % (folder, parameter_folder)):
        os.mkdir("%sTransform_analysis/%s" % (folder, parameter_folder))

    savefolder = "%sTransform_analysis/%s/List_text/" % (folder, parameter_folder)
    if not os.path.isdir(savefolder):
        os.mkdir(savefolder)

    force_list, curve_num_list = PlateauForceList(folder, SD)
    i = 0
    while i < len(force_list):
        Lc = force_list[i]
        if len(Lc) <= 10:
            Lc = zeros(10)
        curve = curve_num_list[i]
        filename = "%sLc_list_%s.txt" %(savefolder, curve)
        saveFile1(filename, Lc, "# Plateau force values (pN) with SD = %g" % SD)
        i = i+1

def locate_tail(tss, force, N = 0, step = 20, SD = 0.017):
    from numpy import std
    '''Locates the final jump.'''
    max_std = 20
    sd = 0
    i = 0
    if N == 0:
        N = len(tss)
    while sd < SD:
        sd = std(force[(N-step-i):(N-i)])
        i = i+1
    return N-i-20

def start_stop(tss, force, SD = 0.017):
    stop = locate_tail(tss, force, SD = SD)
    start = 1
    while force[start] > (min(force)-min(force)*0.01):
        start = start+1
    return start, stop

# More redundant stuff down here

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
    return A, B

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

def curveNumFinder(filename):
    ''' Takes a filename and assuming lineXXXXpontXXXX.txt format, 
        returns the curve number XXXX.XXXX as a string.'''
    import re
    curve =  re.findall('[0-9]{4}', filename)
    curveNum = "%s.%s" % (curve[0], curve[1])
    return curveNum