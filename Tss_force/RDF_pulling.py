''' RDF_pulling.py contains all the functions that load, manipulate, and save a single 
    Ramp, Defl, Force 3-column file as a 2-column file for tss-force.
    Also contains function to save calculate all for single RDF folder.

    Parameters required:
        - name of RDF file 'fileName'
        - location of data set 'folder'
''' 

def curveNumFinder(fileName):
    ''' Takes a filename and assuming XXX.txt format, 
        returns the curve number XXX as a string.
        Aborts and returns error if XXX is not a number.'''
    import re
    curve =  re.search('[a-z]?[0-9]{3}(?=.txt)', fileName)
    curveNum = curve.group(0)
    return curveNum

def loadRDF(fileName):
    ''' Takes a standard SMFS file with ramp, force, and deflection data
        and returns ramp, force, and deflection arrays'''
    from numpy import zeros
    ramp = zeros(6144) # in nm
    defl = zeros(6144) # in nm
    force = zeros(6144) # in pN
    file = open(fileName, 'r')
    file.readline()
    i = 0
    for line in file:
        ramp[i] = float(line.split()[0])
        defl[i] = float(line.split()[1])
        force [i] = float(line.split()[2])
        i = i + 1
    file.close
    return ramp, defl, force

def tss_force(fileName):
    ''' Takes ramp, force and deflection arrays and 
        calculates the aligned and zeroed tss (nm) and force (pN)'''
    from numpy import average
    ramp, defl, force = loadRDF(fileName)
    tss = (ramp-defl)*-1+400
    force = (force-average(force[5837:6144]))*-1
    
    zero=0
    while (force[zero]<0):
        zero=zero+1
    
    tss =tss-tss[zero]
    return [tss, force]

