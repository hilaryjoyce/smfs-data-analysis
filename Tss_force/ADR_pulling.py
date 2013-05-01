def curveNumFinder(filename):
    ''' Takes a filename and assuming lineXXXXpontXXXX.txt format, 
        returns the curve number XXXX.XXXX as a string.'''
    import re
    curve =  re.findall('[0-9]{4}', filename)
    curveNum = "%s.%s" % (curve[0], curve[1])
    return curveNum

def load3Col(fileName):
    '''Takes this weird line break-less data file with 3 columns
    Returns approach, defl, and retract in nanometers'''
    from numpy import zeros
    A = []
    B = []
    C = []
    file = open(fileName, 'r')
    buffer = file.read().split()
    N  = len(buffer)
    A = zeros(N/3)
    B = zeros(N/3)
    C = zeros(N/3)
    i = 0
    k = 0
    while i < len(buffer):
        A[k] = float(buffer[i])
        B[k] = float(buffer[i+1])
        C[k] = float(buffer[i+2])
        k = k + 1
        i = i + 3
    file.close()
    return [A, B, C]

def ADR_to_RDF(fileName, k):
    '''Take a four column file and spring constant and return
    the RDF data in three columns.'''
    from scipy import where
    curve = load3Col(fileName)
    approach = curve[0] # m zsensor in nanometers for approach
    defl = curve[1] # V deflection volts (for approach and retract?)
    retract = curve[2] # m zsensor in nanometers for retract
    force = k*defl
    max_force = where(force == max(force))[0][0]
    #Cut off approach data
    ramp_retract = retract[max_force:-1]*10**9 #nm
    defl_retract = defl[max_force:-1]*10**9 #nm
    force_retract = force[max_force:-1]*10**9 #nN
    return [ramp_retract, defl_retract, force_retract]

def locate_tail(tss, force, N = 0, step = 20, SD = 0.013):
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
        if i == N-25:
            print "reached end"
    return N-i-20

def baseline_fit(tss, force, SD = 0.013):
    from scipy.stats import linregress
    '''Finds a fit to the baseline of the approach curve'''
    #takes approach curve and returns a fit to the last linear regime
    N = locate_tail(tss, force, SD = SD)
    slope, intercept, r, p, se = linregress(tss[N:-1], force[N:-1])
    return slope, intercept

def tss_force(fileName, k, SD = 0.013):
    ''' Adjusts the baseline using the approach slope, and min tss value'''
    RDF = ADR_to_RDF(fileName, k)
    tss = (RDF[0] - RDF[1])*-1
    tss = tss - min(tss)
    force = RDF[2]
    slope, intercept = baseline_fit(tss, force, SD)
    return [tss, force-tss*slope-intercept]
