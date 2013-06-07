''' 
Raw_data_transform.py
Used to transform tss_force data into smoothed or not smoothed curves
by removing adhesion peaks.
'''


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

def all_removeAdhesionPeak(folder):
    '''
    Takes a folder and returns cut_tss_list and cut_force_list,
    which contain tss and force arrays without adhesion peaks.
    '''
    tss_list, force_list = loadTssForce(folder)
    cut_tss_list = []
    cut_force_list = []
    i = 0
    while i < len(tss_list):
        cut_tss, cut_force = remove_adhesion_peak(tss_list[i], force_list[i])
        cut_tss_list.append(cut_tss)
        cut_force_list.append(cut_force)
        i = i+1
    return cut_tss_list, cut_force_list

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

def all_limitTssForce(folder):
    '''
    Gets all the tss and force arrays for a given experiment
    and returns lists of the limited tss and force arrays (removing
    adhesion peaks and baseline).
    '''
    tss_list, force_list = loadTssForce(folder)
    lim_tss_list = []
    lim_force_list = []
    
    for [tss, force] in zip(tss_list, force_list):
        tss_lim, force_lim = limitTssForce(tss, force)
        lim_tss_list.append(tss_lim)
        lim_force_list.append(force_lim)

    return lim_tss_list, lim_force_list

def smoothForce(force, f = 10):
    '''
    Takes a single tss and force array pair and returns 
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