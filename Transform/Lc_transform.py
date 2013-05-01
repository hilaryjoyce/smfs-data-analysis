'''This module deals with take tss-force files and turning them into a list of valid contour lengths.
    Parameters:
        - folder name 'folder'
        - persistence length 'p'
        - minimum valid force 'minForce'
        - minimum valid tss 'minTss'
    Saves into "Transform_analysis" folder (Transform_analysis folders contain any 1D representation of the data!)
'''

def LcList(folder, p = 0.4, minTss = 5, minForce = 45):
    from glob import glob

    tss_force_files = glob('%sTss_Force_text/*.txt' % folder)
    
    Lc_list = []
    curve_num_list = []
    for file in tss_force_files:
        tss, force = load2Col(file)
        tss_lim, force_lim = limitTssForce(tss, force, minTss, minForce)
        Lc_list.append(calcLc(tss_lim, force_lim, p))
        curve_num_list.append(curveNumFinder(file))

    return Lc_list, curve_num_list

def saveLcList(folder, p = 0.4, minTss = 5.0, minForce = 45.0):
    import os

    if not os.path.isdir("%sTransform_analysis/" % folder):
        os.mkdir("%sTransform_analysis/" % folder)

    parameter_folder = "Parameters_p%g_minTss%d_minForce%d/" % (p, int(minTss), int(minForce))
    if not os.path.isdir("%sTransform_analysis/%s" % (folder, parameter_folder)):
        os.mkdir("%sTransform_analysis/%s" % (folder, parameter_folder))

    savefolder = "%sTransform_analysis/%s/List_text/" % (folder, parameter_folder)
    if not os.path.isdir(savefolder):
        os.mkdir(savefolder)

    Lc_list, curve_num_list = LcList(folder, p, minTss, minForce)
    i = 0
    while i < len(Lc_list):
        Lc = Lc_list[i]
        curve = curve_num_list[i]
        filename = "%sLc_list_%s.txt" %(savefolder, curve)
        saveFile1(filename, Lc, "Fitted Contour Lengths (nm) with p = %g" % p)
        i = i+1

def curveNumFinder(fileName):
    ''' Takes a filename and assuming XXX.txt format, 
        returns the curve number XXX as a string.
        Aborts and returns error if XXX is not a number.'''
    import re
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

def limitTssForce(tss, force, minTss=5, minForce=45):
    ''' Calculates tss_lim, force_lim based on input
        parameters of minTss (in nm) and minForce (in pN)'''
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

    force_short = force[(tss > minTss) & (tss < tss_limit)]
    tss_short = tss[(tss>minTss) & (tss < tss_limit)]
    
    force_lim = force_short[(force_short > minForce) & (force_short < 500)]
    tss_lim = tss_short[(force_short > minForce) & (force_short < 500)]
    
    return tss_lim, force_lim

def calcLc(tss_lim, force_lim, p=0.4):  
    ''' Takes the limited tss and force data and an optional 
        persistence length parameter p (nm) (defaulted to 0.4nm) 
        and calculates the contour length array Lc'''
    from numpy import zeros, empty, sqrt
    
    limLength = len(tss_lim)
    Lc = zeros(limLength) # initialize Lc array
    
    alpha = 1.3806488*10**(-23)*298./(p*(10**(-9))) # (KbT/p, units of m)
    
    omega_r = 4*(force_lim*(10**(-12)))/alpha - 3
    omega = empty(limLength, dtype=complex)
    gamma = empty(limLength, dtype=complex)
    Preal = zeros(limLength)

    k=0
    for i in omega_r:
        omega[k] = complex(i,0)
        k=k+1
    k=0
    for w in omega:
        g = (216-w**3+12*sqrt(324-3*w**3))**(1.0/3)
        P = (1-(g+w**2/g - w)/12)**(-1)
        Pr = P.real
        gamma[k] = g
        Preal[k] = Pr
        k=k+1
        
    Lc = Preal*tss_lim
    return Lc

def saveFile1(path, col1, col1_name):
    '''Saves a file with 1 column of data.'''
    file = open(path, 'w')
    file.write("%s\n" % col1_name)
    k  =0
    while k<len(col1):
        el = '%3g\n' % col1[k]
        file.write(el)
        k=k+1
    file.close