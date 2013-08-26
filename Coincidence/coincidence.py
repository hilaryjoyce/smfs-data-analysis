''' Using scipy.signal's fast fourier transform convolution.
    Currently no accounting for shift.
    Parameters:
        - folder specifying the parameters of the inputs (folder)
        - mode of fft (currently just using 'full')
        - shift (currently not even using this parameter)
'''

def absCoincidence(lc, d1, d2, shift='No'):
    '''Runs in one function. Most efficient.'''
    from scipy.signal import fftconvolve
    from numpy import dot, sqrt, where, log10

    dx = lc[1]-lc[0]

    if shift == 0:
        fft_flip = fftconvolve(d1, d2[::-1], mode='valid')
    else:
        fft_flip = fftconvolve(d1, d2[::-1], mode='full')

    if shift == 'No' or shift == 0:
        max_conv = max(fft_flip)
        max_id = where(fft_flip == max_conv)[0][0]
        if shift != 0:        
            max_id = max_id - len(fft_flip)/2
    else:
        shift_id = shift/dx
        shift_left = int(- shift_id + len(fft_flip)/2)
        shift_right = int(shift_id + len(fft_flip)/2)
        max_conv = max(fft_flip[shift_left:shift_right])
        max_id = where(fft_flip == max_conv)[0][0]
        max_id = max_id - len(fft_flip)/2

    max_lc = max_id*dx

    dot_d1 = dot(d1, d1)
    dot_d2 = dot(d2, d2)

    if (dot_d1 <= dot_d2):
        s = 1/dot_d2
    else:
        s = 1/dot_d1

    Gamma = max_conv * s

    return round(Gamma,5), max_id, round(max_lc, int(-1*log10(dx)+2))


def CoincidenceList(folder, shift = 'No', type = 'density'):
    '''
    Calculates a list of the maximum coincidence and corresponding shift 
    for ONE particular shift.
    Folder is the parameter folder starting from "../Data/".
    Type is 'density' or 'force'. 
    '''
    from glob import glob
    if type == 'density':
        files = glob("%sDensity_text/*.txt" % folder)
    else:
        files = glob("%sForce_text/*.txt" % folder)
    eps = 10**(-3) # fixing this since I haven't changed it in so long I forget what it means...

    N = len(files)
    count = 0
    i = 0
    curve1_list = []
    curve2_list = []
    gamma_list = []
    x_max_list = []
    while (i < len(files)):
        curve1 = files[i]
        c1, d1 = load2Col(curve1)
        curveNum1 = curveNumFinder(curve1)
        k = i + 1
        while (k < len(files)):
            curve2 = files[k]
            curveNum2 = curveNumFinder(curve2)
            c2, d2 = load2Col(curve2)
            if (max(d1) == 0 or max(d2) == 0):
                if (max(d1) == 0 and max(d2) == 0):
                    Gamma, u_max, x_max = 1.0, 0, 0
                else:
                    Gamma, u_max, x_max = 0, 0, 0
            else:
                Gamma, u_max, x_max = absCoincidence(c1, d1, d2, shift)
            curve1_list.append(curveNum1)
            curve2_list.append(curveNum2)
            gamma_list.append(Gamma)
            x_max_list.append(x_max)
            k = k+1
        i += 1

    return [curve1_list, curve2_list, gamma_list, x_max_list]

def saveCoincidence(folder, shift = 'No', type='density'):
    '''Tales a folder specifying the PARAMETERS (not Density_text) folder
        and saves a coincidence_report.txt file in a shift folder.sub_tss = where(tss < location + d and tss > location - d)
        Currently this guy only does No shift ('valid' absCoincidence)
        or a full no-shift coincidence for any specified numeric shift.
    '''
    import os
    import sys

    if shift == 'No':
        savefolder = "NoShift/"
    else:
        savefolder = "Shift_%d/" % int(shift)

    if not os.path.isdir(folder+savefolder):
        os.mkdir(folder+savefolder)

    print folder+savefolder

    savefile = "%s%scoincidence_report.txt" % (folder, savefolder)    
    co_matrix = CoincidenceList(folder, shift, type=type)

    file = open(savefile, 'w')
    file.write("# Curve1\tCurve2\tGamma\tBest shift (nm)\n")

    i = 0
    while i < len(co_matrix[0]):
        line = "%s\t%s\t%.3f\t%5.2f\n" % (co_matrix[0][i], co_matrix[1][i], co_matrix[2][i], co_matrix[3][i])
        file.write(line)
        i = i+1

    file.close()

def multipleCoincidence(lc, d1, d2, shift_list):
    '''
    Calculates the coincidence value and associating shift for a list of maximum shifts
    for a single pair of curves.
    '''
    from scipy.signal import fftconvolve
    from numpy import dot, sqrt, where, log10

    dx = lc[1]-lc[0]
    fft_flip = fftconvolve(d1, d2[::-1], mode='full')

    Gamma_list = []
    x_list = []
    u_list = []

    for shift in shift_list:
        if shift == 0:
            max_conv = fft_flip[len(fft_flip)/2]
            max_id = 0
        elif shift == 'No':
            max_conv = max(fft_flip)
            max_id = where(fft_flip == max_conv)[0][0]     
            max_id = max_id - len(fft_flip)/2
        else:
            shift_id = shift/dx
            shift_left = int(- shift_id + len(fft_flip)/2)
            shift_right = int(shift_id + len(fft_flip)/2)
            max_conv = max(fft_flip[shift_left:shift_right])
            max_id = where(fft_flip == max_conv)[0][0]
            max_id = max_id - len(fft_flip)/2

        max_lc = max_id*dx

        dot_d1 = dot(d1, d1)
        dot_d2 = dot(d2, d2)

        if (dot_d1 <= dot_d2):
            s = sqrt(dot_d1/dot_d2)
        else:
            s = sqrt(dot_d2/dot_d1)
        
        Gamma = max_conv / sqrt(dot_d1 * dot_d2) * s
        Gamma_list.append(round(Gamma,5))
        x_list.append(max_id)
        u_list.append(round(max_lc, int(-1*log10(dx)+2)))

    return Gamma_list, x_list, u_list

def MultipleCoincidenceList(folder, type='density', max_x = 600, shift_list = [0, 5, 10, 15, 20, 30, 50, 100, 'No']):
    '''
    Calculates the maximum coincidence and associated shift for a list of maximum shifts
    for an entire folder of curves (param_folder).
    '''
    from glob import glob
    from time import clock
    from sys import exit
    if type == 'density':
        files = glob("%sDensity_text_max%g/*.txt" % (folder, max_x))
    else:
        files = glob("%sForce_text/F*.txt" % folder)
        print files
    eps = 10**(-3)
    N = len(files)
    M = len(shift_list)
    count = 0
    i = 0
    curve1_list = []
    curve2_list = []
    gamma_list = []
    x_max_list = []

    file = open("%sco_run_report.txt" % folder, 'w')
    file.write("First curves completed out of %d:\n" % N)

    while (i < len(files)):
        t1 = clock()
        curve1 = files[i]
        c1, d1 = load2Col(curve1)
        curveNum1 = curveNumFinder(curve1)
        k = i + 1
        while (k < len(files)):
            curve2 = files[k]
            curveNum2 = curveNumFinder(curve2)
            c2, d2 = load2Col(curve2)
            if (max(d1) == 0 or max(d2) == 0):
                u_max, x_max = [0]*M, [0]*M
                if (max(d1) == 0 and max(d2) == 0):
                    Gamma = [1.0]*M
                else:
                    Gamma = [0]*M
            else:
                g = []
                u = []
                x = []
                for shift in shift_list:
                    Gamma, u_max, x_max = multipleCoincidence(c1, d1, d2, shift_list)
            curve1_list.append(curveNum1)
            curve2_list.append(curveNum2)
            gamma_list.append(Gamma)
            x_max_list.append(x_max)
            k = k+1
        i += 1
        t2 = clock()
        print "Completed %d of %d in %g minutes." % (i, N, (t2-t1)/60.0)
        file.write("Completed %d of %d in %g minutes.\n" % (i, N, (t2-t1)/60.0))
    
    file.close()

    return [curve1_list, curve2_list, gamma_list, x_max_list]

def saveMultipleCoincidence(folder, type = 'density', max_x = 600, shift_list = [0, 5, 10, 15, 20, 30, 50, 100, 'No']):
    '''Tales a folder specifying the PARAMETERS (not Density_text) folder
        and saves a coincidence_report.txt file in a shift folder.
    '''
    import os
    import sys

    co_matrix = MultipleCoincidenceList(folder, type, max_x, shift_list)

    if type == 'density':
        co_folder = "%sCoincidence_max%g/" % (folder, max_x)
    else:
        co_folder = "%sCoincidence/" % folder

    if not os.path.isdir(co_folder):
        os.mkdir(co_folder)

    j = 0
    for shift in shift_list:
        if shift == 'No':
            savefolder = "%sNoShift/" % co_folder
        else:
            savefolder = "%sShift_%d/" % (co_folder, int(shift))

        if not os.path.isdir(savefolder):
            os.mkdir(savefolder)

        print savefolder

        savefile = "%scoincidence_report.txt" % savefolder

        file = open(savefile, 'w')
        file.write("# c1\tc2\tGamma\tBest shift (nm)\n")

        i = 0
        N = len(co_matrix[0])
        while i < len(co_matrix[0]):
            line = "%s\t%s\t%.3f\t%5.2f\n" % (co_matrix[0][i], co_matrix[1][i], co_matrix[2][i][j], co_matrix[3][i][j])
            file.write(line)
            i = i+1

        file.close()
        j = j+1

''' 
Function to save multiple coincidence values for different shifts in ONE file without
completely swamping the memory.
'''

def saveAllCoincidence(folder, type='density', max_x = 600, shift_list = [0, 5, 10, 15, 20, 30, 50, 100, 'No']):
    '''
    Calculates the maximum coincidence and associated shift for a list of maximum shifts
    for an entire folder of curves (param_folder).
    '''
    from glob import glob
    from time import clock
    import os
    tstart = clock()
    if type == 'density':
        files = glob("%sDensity_text_max%g/*.txt" % (folder, max_x))
    else:
        files = glob("%sForce_text/*.txt" % folder)

    eps = 10**(-3)
    N = len(files)
    M = len(shift_list)
    count = 0
    i = 0
    
    if type == 'density':
        co_folder = "%sCoincidence_max%g/" %(folder, max_x)
    else:
        co_folder = "%sCoincidence/" % folder
    if not os.path.isdir(co_folder):
        os.mkdir(co_folder)

    print co_folder
    file = open("%sall_coincidence_report.txt" % co_folder, 'w')
    header = "#c1\tc2\t"
    for shift in shift_list:
        header = header + "G%s\ts%s\t" % (str(shift), str(shift))
    header = header + '\n'
    file.write(header)
    
    while (i < len(files)):
        t1 = clock()
        curve1 = files[i]
        c1, d1 = load2Col(curve1)
        curveNum1 = curveNumFinder(curve1)
        k = i + 1
        while (k < len(files)):
            curve2 = files[k]
            curveNum2 = curveNumFinder(curve2)
            c2, d2 = load2Col(curve2)
            if (max(d1) == 0 or max(d2) == 0):
                u_max, x_max = [0]*M, [0]*M
                if (max(d1) == 0 and max(d2) == 0):
                    Gamma = [1.0]*M
                else:
                    Gamma = [0]*M
            else:
                for shift in shift_list:
                    Gamma, u_max, x_max = multipleCoincidence(c1, d1, d2, shift_list)
            j = 0
            line = '%s\t%s\t' % (curveNum1, curveNum2) 
            while j < len(shift_list):
                line = line + '%.3f\t%5.2f\t' % (Gamma[j], x_max[j])
                j = j+1
            line = line + '\n'
            file.write(line)
            k = k+1
        i += 1
        t2 = clock()
        print "Completed %d of %d in %g minutes." % (i, N, (t2-t1)/60.0)    
    
    file.close()
    tfinish = clock()
    
    return 'Completed in %g hours' % ((tfinish-tstart)/3600)

''' 
Separate functions 
------------------
'''


def fft_convolution(lc, d1, d2):
    '''Use fast fourier transform to find max correlation'''
    from scipy.signal import fftconvolve
    from numpy import where
    dx = lc[1]-lc[0]
    fft_flip = fftconvolve(d1, d2[::-1])
    max_conv = max(fft_flip)
    max_id = where(fft_flip == max_conv)[0][0]
    max_id_centered = max_id - len(fft_flip)/2
    max_lc = max_id_centered*dx
    return max_conv, max_id_centered, max_lc

def coincidence(lc, d1, d2):
    '''Gets coincidence value using FFT convolution plugged into
    the cosine correlation function.'''
    from numpy import dot, sqrt
    max_conv, max_id, max_lc = fft_convolution(lc, d1, d2)
    dot_d1 = dot(d1, d1)
    dot_d2 = dot(d2, d2)
    C = max_conv/sqrt(dot_d1*dot_d2)
    return C, max_id, max_lc

def scaling(d1, d2):
    '''Scaling factor'''
    from numpy import dot, sqrt
    dot_d1 = dot(d1, d1)
    dot_d2 = dot(d2, d2)
    if (dot_d1 <= dot_d2):
        s = sqrt(dot_d1/dot_d2)
    else:
        s = sqrt(dot_d2/dot_d1)
    return s

def absCoincidence_step(lc, d1, d2):
    '''Absolute coincidence (cosine correlation * scaling factor) 
    using other functions. Don't use this one for efficiency!'''
    C, max_id, max_lc = coincidence(lc, d1, d2)
    s = scaling(d1, d2)
    Gamma = C * s
    return Gamma, max_id, max_lc

# These are super redundant but OH WELL

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

def split_all_coincidence(param_folder, max_x):
    '''
    Takes an all_coincidence_report.txt IN a Coincidence_max_x folder already
    and returns a folder + coincidence report for each shift.
    '''

    from glob import glob
    import os 
    
    co_folder = '%sCoincidence_max%d/' % (param_folder, max_x)
    co_file = glob('%sall*.txt' % co_folder)[0]

    file = open(co_file, 'r')
    header = file.readline().split()

    # deal with the header
    sg_header = header[2:]
    shift_list = []
    i = 0
    while i < len(sg_header):
        try:
            shift_list.append(int(sg_header[i][1:]))
        except:
            shift_list.append(sg_header[i][1:])
        i = i+2
    print shift_list
    N = len(shift_list)

    shift_data = []
    i = 0
    while i < N:
        shift_data.append([[],[]])
        i = i+1
        
    # break up the rows
    data = []
    curve1_list = []
    curve2_list = []
    for line in file:
        row = line.split()
        curve1_list.append(row[0])
        curve2_list.append(row[1])
        i = 0 
        while i < N:
            shift_data[i][0].append(float(row[i*2+2]))
            shift_data[i][1].append(float(row[i*2+3]))
            i = i+1
    file.close()

    shift_folders = []
    for shift in shift_list:
        if shift == 'No':
            savefolder = "NoShift/"
        else:
            savefolder = "Shift_%d/" % int(shift)
        
        if not os.path.isdir(co_folder+savefolder):
            os.mkdir(co_folder+savefolder)
        
        shift_folders.append(co_folder+savefolder)

    i = 0
    for data in shift_data:
        with open('%scoincidence_report.txt' % shift_folders[i], 'w') as file:
            file.write("# c1\tc2\tGamma\tBest shift (nm)\n")
            k = 0
            while k < len(data[0]):
                line = '%s\t%s\t%.3f\t%5.2f\n' % (curve1_list[k], curve2_list[k], data[0][k], data[1][k])
                file.write(line)
                k = k+1
        i = i+1

    return None


