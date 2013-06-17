''' 
Raw_data_transform.py
Used to transform tss_force data into smoothed or not smoothed curves
by removing adhesion peaks.
'''

def projForceList(folder, dz, pm, max_tss):
    '''
    Takes an experiment folder, loads all raw tss force curves,
    then cuts off the tail and adhesion peak and calculates the 
    projected average onto a master tss array.
    '''
    from glob import glob
    from numpy import arange

    tss_force_files = glob('%sTss_Force_text/*.txt' % folder)

    master_tss = arange(0, max_tss, dz)
    force_list = []
    curve_num_list = []
    
    for file in tss_force_files:
        curve_num_list.append(curveNumFinder(file))
        tss, force = load2Col(file)
        lim_tss, lim_force = limitTssForce(tss, force)
        force_proj = force_projection(lim_tss, lim_force, master_tss, pm)
        pos_force_proj = [max(0, x) for x in force_proj]
        force_list.append(pos_force_proj)

    return force_list, curve_num_list

def saveProjForceList(folder, dz = 0.5, pm = 0.5, max_tss = 400):
    import os
    from numpy import zeros, arange

    if not os.path.isdir("%sTransform_analysis/" % folder):
        os.mkdir("%sTransform_analysis/" % folder)

    parameter_folder = "Parameters_dz%g_pm%g_max%d/" % (dz, pm, max_tss)
    if not os.path.isdir("%sTransform_analysis/%s" % (folder, parameter_folder)):
        os.mkdir("%sTransform_analysis/%s" % (folder, parameter_folder))

    savefolder = "%sTransform_analysis/%s/Force_text/" % (folder, parameter_folder)
    if not os.path.isdir(savefolder):
        os.mkdir(savefolder)

    force_list, curve_num_list = projForceList(folder, dz, pm, max_tss)
    master_tss = arange(0, max_tss, dz)
    i = 0
    while i < len(force_list):
        force = force_list[i]
        if len(force) <= 10:
            force = zeros(10)
        curve = curve_num_list[i]
        filename = "%sForce_list_%s.txt" %(savefolder, curve)
        saveFile2(filename, master_tss, force, 'Tss (nm)', 'projForce (pN)')
        i = i+1

    return "%sTransform_analysis/%s" % (folder, parameter_folder)


'''Functions used overall'''

def limitTssForce(tss, force):
    ''' 
    Calculates tss_lim, force_lim by removing BOTH
    the adhesion peak and the baseline for a single
    tss-force pair of arrays.
    '''
    from numpy import std
    N = len(tss)
    sd=0
    i=1
    while sd < 25:
        sd = std(force[(N-200-i):(N-i)])
        i=i+20
        if (N-200-i) <= 0:
            i = 1
            break

    tss_limit = tss[N-i]
    # LIMIT FORCE AND TSS TO VALID RANGE

    tss_short, force_short = remove_adhesion_peak(tss, force)
    force_lim = force_short[tss_short < tss_limit]
    tss_lim = tss_short[tss_short < tss_limit]

    return tss_lim, force_lim

def av_force(tss, force, location, pm):
    from numpy import average, where, all
    region = where(all(zip(tss < location+pm, tss > location-pm), axis=1))[0]
    if len(region)>0:
        return average(force[region])
    else:
        return 0

def force_projection(tss, force, master_tss, pm):
    from numpy import zeros

    average_force = zeros(len(master_tss))
    i = 0
    while master_tss[i] < tss[0]:
        i = i+1
    while master_tss[i] < tss[-1]:
        average_force[i] = av_force(tss, force, master_tss[i], [pm])
        i = i+1
    return average_force


''' Utility functions '''

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

'''------Separate all functions-----'''
'''Not currently in use'''


def loadTssForce(folder):
    from glob import glob

    tss_force_files = glob('%sTss_Force_text/*.txt' % folder)
    
    tss_list = []
    force_list = []
    for file in tss_force_files:
        tss, force = load2Col(file)
        tss_list.append(tss)
        force_list.append(force)

    return tss_list, force_list

def all_removeAdhesionPeak(tss_list, force_list):
    '''
    Takes a folder and returns cut_tss_list and cut_force_list,
    which contain tss and force arrays without adhesion peaks.
    '''
    cut_tss_list = []
    cut_force_list = []
    i = 0
    while i < len(tss_list):
        cut_tss, cut_force = remove_adhesion_peak(tss_list[i], force_list[i])
        cut_tss_list.append(cut_tss)
        cut_force_list.append(cut_force)
        i = i+1
    return cut_tss_list, cut_force_list

def all_limitTssForce(tss_list, force_list):
    '''
    Gets all the tss and force arrays for a given experiment
    and returns lists of the limited tss and force arrays (removing
    adhesion peaks and baseline).
    '''
    lim_tss_list = []
    lim_force_list = []
    
    for [tss, force] in zip(tss_list, force_list):
        tss_lim, force_lim = limitTssForce(tss, force)
        lim_tss_list.append(tss_lim)
        lim_force_list.append(force_lim)

    return lim_tss_list, lim_force_list

def all_projForce(lim_tss_list, lim_force_list, dz = 0.5, pm = 0.5):
    from numpy import ceil, arange
    
    # Find a good upper limit on tss for data set
    max_tss = 0
    for tss in lim_tss_list:
        if max(tss) > max_tss:
            max_tss = max(tss)
    max_tss = ceil(max_tss/100.0)*100
    
    # Create a master tss array
    master_tss = arange(0, max_tss, dz)
    
    # find the average force arrays
    master_force_list = []
    for [tss, force] in zip(lim_tss_list, lim_force_list):
        master_force_list.append(force_projection(tss, force, master_tss, pm))
    return master_force_list

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

''' ======= Definitions for the variable adhesion peak recognition ======= '''

def sort_by_tss(tss, force):
    from numpy import asarray, argsort
    tss = asarray(tss)
    force = asarray(force)
    args = argsort(tss)
    f = force[args]
    t = tss[args]
    return t, f

def maxes_array(tss, force):
    from numpy import roll, asarray, transpose
    t, f = sort_by_tss(tss, force)
    # roll arrays
    fp = roll(f, 50)
    fn = roll(f, -50)
    # find maximum array
    mat = asarray([fp, f, fn])
    tmat = transpose(mat)
    maxes = [max(row) for row in tmat]
    return maxes, t, f

def original_maxes(tss, force):
    from numpy import zeros
    maxes, t, f = maxes_array(tss, force)
    i = 0
    max_orig = zeros(len(maxes))
    while i < len(maxes):
        if maxes[i] == f[i]:
            max_orig[i] = 1
        i = i+1
    return max_orig, maxes, t, f

def find_last_point(tss, force, d=0.1):
    from numpy import std
    max_orig, maxes, t, f = original_maxes(tss, force)
    # Find last point of adhesion region:
    i = 0
    # skip first 0s
    while i < len(max_orig):
        if max_orig[i] == 1:
            break
        i = i+1
    #print 'zeros', i
    # make sure we're above the zero forces
    while f[i] <= 0:
        i = i+1
        if i == len(f):
            return t, f, 0
    #print i
    # skip the first crazy region (hopefully this doesn't bring us to the end...)
    while std(f[10*i:10*(i+1)]) > 100:
            i = i+1
            if 10*(i+1) == len(f):
                return t, f, 0
    #print 'std', i
    # find region with big tss jumps
    while (i) < len(t):
        if t[i+1]-t[i] > d:
            i = i+1
            break
        i = i+1
    #print 'tss jumps', i
    # find first zeroed region after this
    while any(max_orig[i:i+5] == 1):
        i = i + 1
    #print 'rezeroed', i
    # find first increase in f again
    while f[i+1]-f[i] < 0:
        i = i + 1
    i = i+1
    while f[i+1]-f[i] < 0:
        i = i +1
    #print 'rise', i
    return t, f, i
    
def remove_adhesion_peak(tss, force, d=0.1):
    from numpy import copy
    '''
    Given a (zeroed) force and tss curve, and some delta tss = d,
    Finds the initial adhesion peak region and removes it from the 
    data.
    '''
    # create zeroed force array before last_point
    t, f, i = find_last_point(tss, force, d)
    f_new = copy(f[i:])
    t_new = copy(t[i:])
    return t_new, f_new

''' Gaussian smoothing - not using '''


def smoothForce(force, f = 10):
    '''
    Takes a single tss and force array pair and returns 
    a gaussian filtered curve of the same length.
    '''
    from scipy.ndimage.filters import gaussian_filter1d
    sm_force = gaussian_filter1d(force, f)
    return sm_force

def all_smoothLimForce(folder, f = 10):
    lim_tss_list, lim_force_list = all_limitTssForce(folder)
    smooth_force_list = []
    for force in lim_force_list:
        smooth_force_list.append(smoothForce(force, f=f))
    return smooth_force_list

def all_smoothForce(folder, f = 10):
    tss_list, force_list = all_removeAdhesionPeak(folder)
    smooth_force_list = []
    for force in force_list:
        smooth_force_list.append(smoothForce(force, f=f))
    return smooth_force_list


# Redundant!
# def sortTssForce(folder):
#     ''' 
#     Loads all the tss_force files for a given experiment
#     and returns lists of both sorted by the tss.
#     '''
#     from numpy import argsort
#     tss_list, force_list = loadTssForce(folder)
#     sort_tss_list = []
#     sort_force_list = []
#     for [tss, force] in zip(tss_list, force_list):
#         args = argsort(tss)
#         sort_tss_list.append(tss[args])
#         sort_force_list.append(force[args])
#     return sort_tss_list, sort_force_list