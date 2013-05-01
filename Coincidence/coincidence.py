''' Using scipy.signal's fast fourier transform convolution.
    Currently no accounting for shift.
    Parameters:
        - folder specifying the parameters of the inputs (folder)
        - mode of fft (currently just using 'full')
        - shift (currently not even using this parameter)
'''

def absCoincidence(lc, d1, d2, mode='full'):
    '''Runs in one function. Most efficient.'''
    from scipy.signal import fftconvolve
    from numpy import dot, sqrt, where

    dx = lc[1]-lc[0]
    fft_flip = fftconvolve(d1, d2[::-1], mode='full')
    max_conv = max(fft_flip)
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
    return Gamma, max_id, max_lc

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
