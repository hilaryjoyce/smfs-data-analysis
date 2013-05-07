'''Functions to perform clustering'''

def load_coincidence(fileName):
    #checkign these
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

def get_co_curve_numbers(curve1, curve2):
    '''Takes the curve1 and curve2 columns from the coincidence report
    and returns a list of all the curve numbers.
    Now reduntant because exp_curve_numbers should be identical.'''
    i = 0
    co_curve_numbers = [curve1[0]]
    while curve1[i] == curve1[0]:
      co_curve_numbers.append(curve2[i])
      i = i+1
    return co_curve_numbers

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
    return root, node_list

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
    flat cluster ids (M).'''
    from scipy.cluster.hierarchy import leaders
    T = flat_cluster(Z, cut_co)
    node_list = hierarchical_tree(Z)
    (L, M) = leaders(Z, T)
    return L, M

'''Functions for organizing flat cluster information'''

def flat_sizes(flat_cluster):
    '''Given a flat cluster, return a dictionary of 
    cluster ids and their sizes.'''
    N = max(flat_cluster)
    ids = range(1,N+1)
    flat_cluster = list(flat_cluster)
    count = [flat_cluster.count(x) for x in ids]
    return ids, count

def flat_ordered_by_size(flat_cluster):
    '''Given a flat cluster, return a dictionary of 
    cluster ids and their sizes, ordered from largest to smallest
    cluster..'''
    from numpy import asarray, argsort
    ids, count = flat_sizes(flat_cluster)
    ids = asarray(ids)
    count = asarray(count)
    args = argsort(count)
    return ids[args][::-1], count[args][::-1]

def print_sizes(flat_cluster, minSize=1):
    '''Print cluster ids and sizes from largest to the
    minimum specified cluster size (default = 1)'''
    ids, count = flat_ordered_by_size(flat_cluster)
    i = 0
    print "id\tsize"
    while (i < len(count)) and (count[i] >= minSize):
        print "%d\t%d" % (ids[i], count[i])
        i = i+1

'''Little internal guys'''

def cluster_ids(flat_cluster):
    return range(1, max(flat_cluster)+1)

def _map_pair_index(flat_cluster, c1, c2):
    '''Takes the indexes of two curves and the index of each curve
    returns the pair index in the coincidence report.'''
    N = len(flat_cluster)
    if c1 > c2:
        c2, c1 = c1, c2
    pair_index = sum(range(N)[::-1][0:c1]) + c2 - c1 - 1
    return pair_index

def curves_in_cluster(flat_cluster, id):
    '''Given a flat cluster and a cluster id, returns
    a list of all the curve indexes in that cluster'''
    index = range(0, len(flat_cluster))
    return [x for x in index if flat_cluster[x]==id]

def curves_in_clusters(flat_cluster):
    '''Given a flat cluster, returns a dictionary of each
    cluster id and the indexes of all curves in that cluster.'''
    ids = cluster_ids(flat_cluster)
    cluster_curves = [curves_in_cluster(flat_cluster, i) for i in ids]
    return ids, cluster_curves

def map_curve_to_report(flat_cluster, c1, c2, report_list):
    '''Takes a curve number, a pair of curve indexes, and a list
    of the length of the coincidence report and returns
    the value from that list for that curve pair.'''
    co_curve_indexes = get_co_curve_numbers(flat_cluster)
    return report_list[_map_pair_index(co_curve_indexes, c1, c2)]

def map_curve_to_property(curve_index, property_list):
    '''Takes a curve number (indexed) and an ordered property
    list and returns the property for that curve.'''
    return property_list[curve]

'''Functions for looking at individual clusters'''
def map_cluster_to_curves(cluster_number, flat_cluster):
    '''Takes a cluster number and a flattened cluster and
    returns a list of the indexes in that cluster.'''
    index = range(0, len(flat_cluster))
    return [x for x in index if flat_cluster[x]==cluster_number]

def map_cluster_to_pairs(cluster_number, flat_cluster):
    '''Takes a cluster number, a flattened cluster, and the index of 
    all curves and returns a list of the pair indices in the co report.'''
    index = range(0,len(flat_cluster))
    return _map_pair_index_list(index, map_cluster_to_curves(cluster_number, flat_cluster))

def cluster_curve_properties(cluster_number, flat_cluster, prop_list):
    '''Takes a cluster number, a flattened cluster, the curve index, and a list
    (of length number of curves) with some property of each curve in it, and returns 
    a list of only '''
    return [prop_list[i] for i in map_cluster_to_curves(cluster_number, flat_cluster)]    

def cluster_pair_report(cluster_number, flat_cluster, rep_list):
    '''Takes a cluster number, a flattened cluster, 
    and any coincidence report array and returns a list of those values.'''
    index = range(0,len(flat_cluster))
    cluster_curve_indexes = map_cluster_to_curves(cluster_number, flat_cluster)
    return [map_curve_to_report(index, c1, c2, rep_list) \
    for c1 in cluster_curve_indexes for c2 in cluster_curve_indexes if c2 > c1]

'''Shiting and average behaviour'''

'''Actual functions I will use to manipulate clusters for investigation.'''

def map_cluster_to_files(parameter_folder, flat_cluster, id, file_type = 'None'):
    curve_ids = curves_in_cluster(flat_cluster, id)
    files = leaf_files(parameter_folder, file_type)
    return [files[i] for i in curve_ids]

def cluster_density_curves(parameter_folder, flat_cluster, id):
    files = map_cluster_to_files(parameter_folder, flat_cluster, id, file_type = 'density')
    lc_density_list = []
    for file in files:
        lc_density_list.append(load2Col(file))
    return lc_density_list

def cluster_shift_values(flat_cluster, id, shift_array):
    shift_list = cluster_pair_report(id, flat_cluster, shift_array)
    return shift_list[0:len(flat_cluster)]

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
    return [A, B]








