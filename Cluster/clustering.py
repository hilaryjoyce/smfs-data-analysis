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

def leaf_files(param_folder, file_type = 'None'):
    import re
    from glob import glob
    if file_type == 'tss_force':
        folder = re.search('\S+(?=Transform_analysis)', param_folder).group(0)
        tss_force_folder = folder+'Tss_Force_text/'
        return glob(tss_force_folder+"*.txt")
    elif file_type == 'density':
        density_folder = param_folder + "Density_text/"
        return glob(density_folder+"*.txt")
    elif file_type == 'list':
        list_folder = param_folder + "List_text/"
        return glob(list_folder+"*.txt")
    elif file_type == 'RDF':
        folder = re.search('\S+(?=Transform_analysis)', param_folder).group(0)
        RDF_folder = folder + "RDF_text/"
        return glob(RDF_folder+"*.txt")
    else:
        return "Specify file_type = {'tss_force', 'density', 'list'}"

def hierarchical(co_report_file):
    '''Takes a coincidence_report text file location and returns
    the scipy.cluster object.'''
    from scipy.cluster.hierarchy import linkage
    curve1, curve2, co, shift = load_coincidence(co_report_file)
    Z = linkage(1-co, method='complete')
    return Z

def hierarchical_tree(Z):
    '''Takes a linkage matrix and returns the pre-order traversal list of 
    all the ClusterNodes.'''
    from scipy.cluster.hierarchy import to_tree
    root, node_list = to_tree(Z, rd = True)
    return node_list

def flat_cluster(Z, cut_co):
    '''Takes a full linkage matrix Z and returns the flattened cluster
    array T with minimum coincidence of cut_co'''
    from scipy.cluster.hierarchy import fcluster
    t = 1 - cut_co
    T = fcluster(Z, t, criterion = 'distance')
    return T

def flat_cluster_roots(Z, cut_co):
    '''Takes a linkage array and returns a list of the ClusterNodes
    that form the nodes for that flat cluster (N) and their corresponding 
    flat cluster ids (M)'''
    from scipy.cluster.hierarchy import leaders
    T = flat_cluster(Z, cut_co)
    node_list = hierarchical_tree(Z)
    L, M = leaders(Z, T)
    return L, M


