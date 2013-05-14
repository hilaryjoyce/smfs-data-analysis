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
    ''' 
    An CoAnalysis class for representing all the information about all
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

    def __str__(self):
        l1 = self.parameter_folder + '\n'
        l2 = 'Hierarchical clustering at max_shift = %s\n' % str(self.max_shift)
        l3 = 'Number of curves = %d\n' % len(self.list_curve_names())
        return l1+l2+l3

    def get_coincidence_array(self):
        return self.co_array

    def get_shift_array(self):
        return self.shift_array

    def get_curve1_array(self):
        return self.curve1

    def get_curve2_array(self):
        return self.curve2

    def list_density_files(self):
        return self.density_files

    def list_tss_force_files(self):
        return self.tss_force_files

    def list_curve_names(self):
        '''
        Returns a list of all the string names of the files used in
        this clustering. 
        '''
        return [curveNumFinder(file) for file in self.tss_force_files]

    def list_curve_indexes(self):
        '''
        Returns a list [0,...,N-1] where N is the number of curves used in
        this clustering.
        '''
        return range(len(self.tss_force_files))

    def hierarchical_cluster(self):
        '''
        Returns the hierarchical linkage array Z, clustered
        using the complete clustering method in scipy.cluster.hierarchy
        '''
        from scipy.cluster.hierarchy import linkage
        Z = linkage(1-self.co_array, method='complete')
        return Z

    def plot_dendrogram(self):
        '''
        Plots the dendragram visualization of Z and returns
        the dendrogram object 'dendro'.
        '''
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

    # Adding some functions to return pair information for any two curves

    def _map_pair_to_report(self, c1, c2, mode):
        ''' 
        Takes two curves (by index 0:N) and returns the corresponding
        index in the pair report specified by mode ('co', 'shift', 'curve1', 'curve2').
        '''
        curve_indexes = self.list_curve_indexes()
        if c1 > c2:
            c2, c1 = c1, c2
        pair_index = sum(curve_indexes[::-1][0:c1]) + c2 - c1 - 1
        if mode == 'curve1':
            return self.get_curve1_array()[pair_index]
        if mode == 'curve2':
            return self.get_curve2_array()[pair_index]
        if mode == 'co':
            return self.get_coincidence_array()[pair_index]
        if mode == 'shift':
            return self.get_shift_array()[pair_index]
        else:
            print 'Invalid mode, try {co, shift, curve1, curve2}'
            return None

    def map_pair_to_shift(self, c1, c2):
        '''Given a pair of curves by index number, returns the best shift.'''
        return self._map_pair_to_report(c1, c2, mode='shift')

    def map_pair_to_coincidence(self, c1, c2):
        '''Given a pair of curves by index number, returns the best coincidence.'''
        return self._map_pair_to_report(c1, c2, mode='co')

    def map_pair_to_indexes(self, c1, c2):
        '''Given a pair of curves by index number, returns a tuple of the curve names.'''
        if c1 < c2:
            return self._map_pair_to_report(c1, c2, mode='curve1'), self._map_pair_to_report(c1, c2, mode='curve2')
        else:
             return self._map_pair_to_report(c1, c2, mode='curve2'), self._map_pair_to_report(c1, c2, mode='curve1')

class FlatClusters(CoAnalysis):

    def __init__(self, parameter_folder, max_shift, co_cut):
        super(FlatClusters, self).__init__(parameter_folder, max_shift)
        self.co_cut = co_cut

        Z = self.hierarchical_cluster()
        self.T = flatten(Z, 1-self.co_cut)
        self.number = max(self.T)

    def get_parameter_folder(self):
        return self.parameter_folder

    def get_min_coincidence(self):
        '''
        Returns the minimum coincidence value of this flattening.
        '''
        return self.co_cut

    def list_flat_cluster(self):
        '''
        Returns the list representing the cluster id of each curve 
        at this flattening.
        '''
        return list(self.T)

    def get_number_of_clusters(self):
        '''
        Returns the number of clusters at this flattening.
        '''
        return self.number

    def list_cluster_ids(self):
        ''' 
        Retruns a list of [1, ..., N] where N is the number of clusters
        representing the ids of each cluster in the flattened set of clusters.
        '''
        return range(1,self.number+1)

    def list_cluster_sizes(self):
        ''' 
        Returns a pre-ordered list of the sizes of each cluster.
        Cluster n has size in this list at position n-1.
        '''
        return [list(self.T).count(x) for x in self.list_cluster_ids()]

    def ordered_clusters(self):
        '''
        Returns a tuple of the cluster ids array ordered by size, 
        and the correspondingly ordered array of the curve sizes.
        '''
        from numpy import asarray, argsort
        sizes = asarray(self.list_cluster_sizes())
        clusters = asarray(self.list_cluster_ids())
        args = argsort(sizes)[::-1]
        return clusters[args], sizes[args]

    def largest_clusters(self, minSize = 1):
        '''
        Returns a tuple of the cluster ids array ordered by size
        and the correspondingly ordered array of the curve sizes,
        for every cluster at least as large as the minSize.
        '''
        clusters, sizes = self.ordered_clusters()
        i = 0
        largest_sizes = []
        largest_clusters = []
        while i < len(sizes) and sizes[i] >= minSize:
            largest_sizes.append(sizes[i])
            largest_clusters.append(clusters[i])
            i = i+1
        return largest_clusters, largest_sizes

    def curves_by_cluster(self):
        '''
        Returns a list of the lists of all the curves (by index) in
        each cluster (cluster n at postion n-1).
        '''
        flat = self.list_flat_cluster()
        ids = self.list_curve_indexes()
        cluster_list = self.list_cluster_ids()
        clusters = []
        for c in cluster_list:
            clusters.append([id for id in ids if flat[id] == c])
        return clusters

    def __str__(self):
        l1 = self.parameter_folder + '\n'
        l2 = 'Flatted cluster at minimum coincidence %g\n' % self.get_min_coincidence()
        l3 = '%d clusters of %d curves\n' % (self.get_number_of_clusters(), len(self.list_flat_cluster()))
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

        self.min_coincidence = flat_cluster.get_min_coincidence()

    def __str__(self):
        flat = self.flat
        l1 = flat.get_parameter_folder() + '\n'
        l2 = 'Flatted cluster at minimum coincidence %g\n' % self.get_min_coincidence()
        l3 = "Cluster %d with %d curves." %(self.cluster_number, len(self.indexes))
        return l1+l2+l3

    def list_curve_indexes(self):
        '''
        Returns a list of all the curve indices in this cluster.
        '''
        return self.indexes

    def get_flat_cluster(self):
        '''
        Returns the flat clustering object associated with this cluster.
        '''
        return self.flat

    def get_min_coincidence(self):
        '''
        Returns the minimum coincidence of the curves in this cluster.
        '''
        return self.min_coincidence

    def get_cluster_size(self):
        '''
        Returns the size of the cluster.
        '''
        return len(self.indexes)

    def list_curve_names(self):
        '''
        Return the string names for each curve in this cluster.
        '''
        flat_names = self.flat.list_curve_names()
        return [flat_names[i] for i in self.indexes]

    def list_density_files(self):
        '''
        Returns a list of the locations of the density files for this cluster.
        '''
        flat_files = self.flat.list_density_files()
        return [flat_files[i] for i in self.indexes]

    def list_tss_force_files(self):
        '''
        Returns a list of the locations of the tss_force files for this cluster.
        '''
        flat_files = self.flat.list_tss_force_files()
        return [flat_files[i] for i in self.indexes]

    def list_cluster_coincidences(self):
        flat = self.get_flat_cluster()
        indexes = self.list_curve_indexes()
        coincidence = []
        i = 0
        while i < len(indexes):
            k = i+1
            while k < len(indexes):
                coincidence.append(flat.map_pair_to_coincidence(indexes[i], indexes[k]))
                k = k+1
            i = i+1
        return coincidence

    def list_cluster_shift(self):
        flat = self.get_flat_cluster()
        indexes = self.list_curve_indexes()
        shift = []
        i = 0
        while i < len(indexes):
            k = i+1
            while k < len(indexes):
                shift.append(round(flat.map_pair_to_shift(indexes[i], indexes[k]),4))
                k = k+1
            i = i+1
        return shift

    def get_Lc_density_arrays(self):
        density_files = self.list_density_files()
        Lc_density_list = []
        for file in density_files:
            Lc_density_list.append(load2Col(file))
        return Lc_density_list

    def get_shifted_Lc_density_arrays(self):
        Lc_density_list = self.get_Lc_density_arrays()
        shift_list = self.list_cluster_shift()

