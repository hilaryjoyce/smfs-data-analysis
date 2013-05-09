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

def curve_files(param_folder, file_type = 'None'):
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
    Z = hier.linkage(1-co, method='complete')
    return Z

def flatten(Z, t):
    from scipy.cluster.hierarchy import fcluster
    return fcluster(Z, t, criterion='distance')

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

''' --------------------- '''

class CoAnalysis(object):
    ''' An CoAnalysis class for representing all the information about all
    the curves collected for a given experiment and set of analysis parameters. 

    Parameters
    ----------

    parameter_folder :  the full string reference to the location of the parameter folder in
                        the Transform_analysis folder that contains this set of data.

    max_shift :         the maximum shift value used for this coincidence matrix
    '''
    def __init__(self, parameter_folder, max_shift):
        self.parameter_folder = parameter_folder
        self.max_shift = max_shift

        if max_shift == 'No':
            shift_folder = "NoShift/"
        else:
            shift_folder = "Shift_%d/" % max_shift

        # File names for associated parameters
        self.coincidence_file = "%s%scoincidence_report.txt" % (parameter_folder, shift_folder)
        self.tss_force_files = curve_files(parameter_folder, file_type = 'tss_force')
        self.density_files = curve_files(parameter_folder, file_type = 'density')

        # Load coincidence file
        curve1, curve2, co_array, shift_array = load_coincidence(self.coincidence_file)
        self.co_array = co_array
        self.shift_array = shift_array
        self.curve1 = curve1
        self.curve2 = curve2

    def coincidence_array(self):
        return self.co_array

    def shift_array(self):
        return self.shift_array

    def curve_names(self):
        return [curveNumFinder(file) for file in self.tss_force_files]

    def curve_indexes(self):
        return range(len(self.tss_force_files))

    def hierarchical_cluster(self):
        from scipy.cluster.hierarchy import linkage
        Z = linkage(1-self.co_array, method='complete')
        return Z

    def plot_dendrogram(self):
        from scipy.cluster.hierarchy import dendrogram
        import matplotlib.pyplot as plt

        h_clustering = self.hierarchical_cluster()

        dendro = dendrogram(h_clustering, no_labels=True, count_sort=False, orientation="left");
        plt.title("Clustering Diagram for N = 3201", fontsize = 14)
        plt.xlabel("Coincidence Metric ($\Gamma$)", fontsize = 14)
        plt.ylabel("Clusters", fontsize = 14)
        plt.xticks([1, 0.8, 0.6, 0.4, 0.2, 0], [0, 0.2, 0.4, 0.6, 0.8, 1])
        plt.xticks(fontsize=14)
        return dendro

class FlatClusters(CoAnalysis):

    def __init__(self, parameter_folder, max_shift, co_cut):
        super(FlatClusters, self).__init__(parameter_folder, max_shift)
        self.co_cut = co_cut

        Z = self.hierarchical_cluster()
        self.T = flatten(Z, 1-self.co_cut)
        self.number = max(self.T)

    def minimum_coincidence(self):
        return self.co_cut

    def flat_cluster(self):
        return self.T

    def number_of_clusters(self):
        return self.number

    def cluster_list(self):
        return range(1,self.number+1)

    def cluster_sizes(self):
        return [list(self.T).count(x) for x in self.cluster_list()]

    def largest_clusters(self, minSize = 1):
        clusters, sizes = self.ordered_clusters()
        i = 0
        largest_sizes = []
        largest_clusters = []
        while i < len(sizes) and sizes[i] >= minSize:
            largest_sizes.append(sizes[i])
            largest_clusters.append(clusters[i])
            i = i+1
        return largest_clusters, largest_sizes

    def ordered_clusters(self):
        from numpy import asarray, argsort
        sizes = asarray(self.cluster_sizes())
        clusters = asarray(self.cluster_list())
        args = argsort(sizes)[::-1]
        return clusters[args], sizes[args]

    def curves_by_cluster(self):
        flat = self.flat_cluster()
        ids = self.curve_indexes()
        cluster_list = self.cluster_list()
        clusters = []
        for c in cluster_list:
            clusters.append([id for id in ids if flat[id] == c])
        return clusters

    def __str__(self):
        l1 = self.parameter_folder + '\n'
        l2 = 'Flatted cluster at minimum coincidence %g\n' % self.minimum_coincidence()
        l3 = '%d clusters of %d curves\n' % (self.number_of_clusters(), len(self.flat_cluster()))
        return l1+l2+l3

class Cluster:

    '''
        Given a given flat clustering and a number for a cluster in that flat clustering,
        returns a single 'Cluster' object that contains information about that cluster
        and its curves.

        Parameters
        ---------
        
        flat_cluster : a flattened clustering object (FlatClusters) from a hierarchical cluster object
        
        cluster_number : a number from 1 to N where N is the number of clusters in flat_cluster

    ''' 

    def __init__(self, flat_cluster, cluster_number):
        
        self.cluster_number = cluster_number
        self.flat = flat_cluster
        self.indexes = flat_cluster.curves_by_cluster()[cluster_number-1]

        self.min_coincidence = flat_cluster.minimum_coincidence()

    def list_curve_indexes(self):
        return self.indexes

    def get_min_coincidence(self):
        return self.min_coincidence

    def get_cluster_size(self):
        return len(self.indexes)

    def list_curve_names(self):
        '''Return the string names for each curve in this cluster.'''
        flat_names = self.flat.curve_names()
        curve_names = []
        curve_indexes = self.curve_indexes()
        i = 0
        while i < len(curve_indexes):
            curve_names.append(flat_names[curve_indexes[i]])
            i = i+1
        return curve_names


