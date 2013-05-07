'''Functions to perform clustering'''

def load_coincidence(fileName):
    '''Loads a single coincidence_report.txt file.
    Returns curve1, curve2, coincidence, shift as lists.'''
    from numpy import asarray
    import re
    curve1 = []
    curve2 = []
    co = []
    shift = []
    file = open(fileName, 'r')
    for line in file:
        try:
            co.append(float(line.split()[2]))
            shift.append(float(line.split()[3]))
            c1 = line.split()[0]
            c2 = line.split()[1]
            c1 = re.search('([a-z]?)\d+', c1).group(0)
            c2 = re.search('([a-z]?)\d+', c2).group(0)
            curve1.append(c1)
            curve2.append(c2)
        except:
            continue
    file.close
    co = asarray(co)
    shift = asarray(shift)
    return curve1, curve2, co, shift

def hierarchical(co_report_file):
    '''Takes a coincidence_report text file location and returns
    the scipy.cluster object.'''
    from scipy.cluster.hierarchy import linkage
    curve1, curve2, co, shift = load_coincidence(co_report_file)
    Z = hier.linkage(1-co, method='complete')
    return Z

def flatten(Z, t):
    from scipy.cluster.hierarchy import fcluster
    return fcluster(Z, t, criterion='distance')

class FlatCluster:

    def __init__(self, Z, t):
        
