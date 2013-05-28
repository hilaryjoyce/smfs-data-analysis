'''Functions to perform clustering'''

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

        # Find correct folder for coincidence data
        if max_shift == 'No':
            shift_folder = "NoShift/"
        else:
            shift_folder = "Shift_%d/" % max_shift

        self.shift_folder = shift_folder

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

        # Let's caculate Z on initiation, since we will always need it!
        self.Z = self.hierarchical_cluster()

    def __str__(self):
        l1 = self.parameter_folder + '\n'
        l2 = 'Hierarchical clustering at max_shift = %s\n' % str(self.max_shift)
        l3 = 'Number of curves = %d\n' % len(self.list_curve_names())
        return l1+l2+l3

    def get_sample_size(self):
        '''Returns the number of initial data points.'''
        return len(self.tss_force_files)

    def get_coincidence_array(self):
        '''Returns the array of coincidence values.'''
        return self.co_array

    def get_shift_array(self):
        '''Returns the array of shift values'''
        return self.shift_array

    def get_curve1_array(self):
        '''Returns the array of curve1 names.'''
        return self.curve1

    def get_curve2_array(self):
        '''Returns the array of curve2 names.'''
        return self.curve2

    def get_parameter_folder(self):
        '''Returns the string of the parameter folder.'''
        return self.parameter_folder

    def get_hierarchical_cluster(self):
        '''Returns the flattened clustering matrix Z.'''
        return self.Z

    def list_density_files(self):
        '''Returns a list of all the density files.'''
        return self.density_files

    def list_tss_force_files(self):
        '''Returns a list of all the tss_force files.'''
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
        from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
        import matplotlib.pyplot as plt

        cpool = ["#1F78B4", "#E31A1C", "#A6CEE3", "#FB9A99", "#7BCCC4", "#B2DF8A", "#33A02C", "#02818A", "#FF7F00", "#FDBF6F", "#CAB2D6", "#6A3D9A", "#BFD3E6", "#8C96C6"]
        set_link_color_palette(cpool)
        h_clustering = self.Z

        dendro = dendrogram(h_clustering, no_labels=True, 
            count_sort=True, orientation="left");
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

class FlatClusters:

    def __init__(self, co_analysis, co_cut):
        ''' 
        Instead of inheriting and thus recalculating a CoAnalysis object
        we will just require one.
        '''    
        self.co_analysis = co_analysis
        self.co_cut = co_cut
        self.parameter_folder = co_analysis.get_parameter_folder()

        Z = co_analysis.get_hierarchical_cluster()
        self.T = flatten(Z, 1-self.co_cut)
        self.number = max(self.T)

    def __str__(self):
        l1 = self.parameter_folder + '\n'
        l2 = 'Flatted cluster at minimum coincidence %g\n' % self.get_min_coincidence()
        l3 = '%d clusters of %d curves\n' % (self.get_number_of_clusters(), len(self.list_flat_cluster()))
        return l1+l2+l3
    
    def get_co_analysis(self):
        ''' 
        Returns the CoAnalysis object used to perform this flattening.
        '''
        return self.co_analysis

    def get_parameter_folder(self):
        '''
        Returns the parameter folder for this Transform analysis.
        '''
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

    def print_clusters_by_size(self, minSize = 1):
        largest_clusters, largest_sizes = self.largest_clusters(minSize=minSize)
        i = 0
        print "Cluster\tSize"
        while i < len(largest_clusters):
            print "%d\t%d" % (largest_clusters[i], largest_sizes[i])
            i = i+1

    def curves_by_cluster(self):
        '''
        Returns a list of the lists of all the curves (by index) in
        each cluster (cluster n at postion n-1).
        '''
        flat = self.list_flat_cluster()
        coa = self.co_analysis
        ids = coa.list_curve_indexes()
        cluster_list = self.list_cluster_ids()
        clusters = []
        for c in cluster_list:
            clusters.append([id for id in ids if flat[id] == c])
        return clusters

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
        self.co_analysis = flat_cluster.get_co_analysis()
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

    def get_cluster_size(self):
        return len(self.indexes)

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
        coa = self.co_analysis
        flat_names = coa.list_curve_names()
        return [flat_names[i] for i in self.indexes]

    def list_density_files(self):
        '''
        Returns a list of the locations of the density files for this cluster.
        '''
        coa = self.co_analysis
        flat_files = coa.list_density_files()
        return [flat_files[i] for i in self.indexes]

    def list_tss_force_files(self):
        '''
        Returns a list of the locations of the tss_force files for this cluster.
        '''
        coa = self.co_analysis
        flat_files = coa.list_tss_force_files()
        return [flat_files[i] for i in self.indexes]

    def list_cluster_coincidences(self):
        coa = self.co_analysis
        indexes = self.list_curve_indexes()
        coincidence = []
        i = 0
        while i < len(indexes):
            k = i+1
            while k < len(indexes):
                coincidence.append(round(coa.map_pair_to_coincidence(indexes[i], indexes[k]),5))
                k = k+1
            i = i+1
        return coincidence

    def list_cluster_shifts(self):
        coa = self.co_analysis
        indexes = self.list_curve_indexes()
        shift = []
        i = 0
        while i < len(indexes):
            k = i+1
            while k < len(indexes):
                shift.append(round(coa.map_pair_to_shift(indexes[i], indexes[k]),4))
                k = k+1
            i = i+1
        return shift

    def list_initial_shifts(self):
        shift_list = self.list_cluster_shifts()
        initial_shifts = [0] + shift_list[0:self.get_cluster_size()-1]
        return initial_shifts

    def get_Lc_density_arrays(self):
        density_files = self.list_density_files()
        Lc_density_list = []
        for file in density_files:
            Lc_density_list.append(load2Col(file))
        return Lc_density_list

    def get_shifted_Lc_density_arrays(self):
        Lc_density_list = self.get_Lc_density_arrays()
        shift_list = self.list_cluster_shifts()

    def plot_cluster(self, alpha = 0.5, max_x = 150):
        import matplotlib.pyplot as plt
        from numpy import average, arange
        Lc_density_list = self.get_Lc_density_arrays()
        initial_shifts = self.list_initial_shifts()
        av_shift = average(initial_shifts)
        co = self.get_min_coincidence()

        i = 0
        for Lc_density in Lc_density_list:
            plt.plot(Lc_density[0] + initial_shifts[i] - av_shift, Lc_density[1], 'k-', alpha=alpha)
            i = i+1
        plt.xlim(0,max_x)
        plt.yticks([])
        max_y = plt.axis()[3]
        #plt.text(max_x/1.5, max_y/1.5, '$\Gamma\geq$ %.3f' % co, size = 16)
        #plt.xlabel("Contour Length (nm)")
        #plt.ylabel("Density")
        #plt.title("Cluster of %d curves at co = %g" % (self.get_cluster_size(), self.get_min_coincidence()))
        plt.text(max_x/1.5, max_y/1.5, '%d curves' % self.get_cluster_size(), size = 12)
        return plt.gcf()

class SubCluster:

    def __init__(self, flat_cluster, cluster = None, curve_list = None):
        from scipy.cluster.hierarchy import to_tree
        from numpy import asarray, sort

        self.flat = flat_cluster # FlatClusters object
        self.co_analysis = self.flat.get_co_analysis() #CoAnalysis object
        self.cluster = cluster #Cluster object

        if not cluster == None:
            self.curve_list = cluster.list_curve_indexes()
        else:
            self.curve_list = curve_list

        self.Z = self.co_analysis.get_hierarchical_cluster()

        root = to_tree(self.Z) # root of entire cluster!
        curves = asarray(self.curve_list) # list of curves in this cluster

        # Get the cluster node that corresponds to the curves in the cluster above
        self.cluster_node = get_cluster_node(root, root.left, root.right, curves)
        self.id = self.cluster_node.get_id()

        # Get the right and left cluster nodes
        self.left = self.cluster_node.left
        self.right = self.cluster_node.right

        # Get the left and right cluster lists
        self.left_list = sort(any_pre_order(root, self.left))
        self.right_list = sort(any_pre_order(root, self.right))

    def subcluster_properties(self, property = None):
        ''' 
        Returns a tuple of the left, right subcluster properties asked for in "property".
        property = {size, indexes, names, density_files, Lc_density, shifts, initial_shifts}
        '''
        if property == 'size':
            return (self.get_cluster_size('l'), self.get_cluster_size('r'))
        elif property == 'indexes':
            return (self.list_curve_indexes('l'), self.list_curve_indexes('r'))
        elif property == 'names':
            return (self.list_curve_names('l'), self.list_curve_names('r'))
        elif property == 'density_files':
            return (self.list_density_files('l'), self.list_density_files('r'))
        elif property == 'Lc_density':
            return (self.get_Lc_density_arrays('l'), self.get_Lc_density_arrays('r'))
        elif property == 'shifts':
            return (self.list_curve_shifts('l'), self.list_curve_shifts('r'))
        elif property == "initial_shifts":
            return (self.list_initial_shifts('l'), self.list_initial_shifts('r'))
        else:
            return ('left', 'right')

    def get_cluster_node(self, subcluster = 'subcluster'):
        ''' 
        Returns the ClusterNode object referring to either the supercluster
        or the left and right subclusters.
        ''' 
        if subcluster == 'left' or subcluster == 'l':
            return self.left
        elif subcluster == 'right' or subcluster == 'r':
            return self.right
        else:
            return self.cluster_node

    def get_cluster_size(self, subcluster):
        '''Return size of the left or right subclusters'''
        if subcluster == 'left' or subcluster == 'l':
            indexes = self.left_list
        elif subcluster == 'right' or subcluster == 'r':
            indexes = self.right_list
        return len(indexes)

    def get_cluster_coincidence(self, subcluster):
        '''
        Returns the minimum coincidence for the left or right clusters.
        '''
        if subcluster == 'left' or subcluster == 'l':
            cluster = self.left
        if subcluster == 'right' or subcluster == 'r':
            cluster = self.right
        return 1-cluster.dist

    def list_curve_indexes(self, subcluster):
        '''
        Returns an ordered list of curve ids.
        '''
        if subcluster == 'left' or subcluster == 'l':
            indexes = self.left_list
        if subcluster == 'right' or subcluster == 'r':
            indexes = self.right_list
        return indexes

    def list_curve_names(self, subcluster):
        '''
        Returns a list of the cluster names in the left or right clusters.
        '''
        if subcluster == 'left' or subcluster == 'l':
            indexes = self.left_list
        if subcluster == 'right' or subcluster == 'r':
            indexes = self.right_list
        coa = self.co_analysis
        names = coa.list_curve_names()
        return [names[i] for i in indexes]

    def list_density_files(self, subcluster):
        '''
        Returns a list of the density files in the specified subcluster.
        '''
        if subcluster == 'left' or subcluster == 'l':
            indexes = self.left_list
        if subcluster == 'right' or subcluster == 'r':
            indexes = self.right_list
        coa = self.co_analysis
        flat_files = coa.list_density_files()
        return [flat_files[i] for i in indexes]

    def get_Lc_density_arrays(self, subcluster):
        '''
        Returns a list of pairs of the Lc_arrays and density_arrays for 
        the specified subcluster {left, right}.
        '''
        density_files = self.list_density_files(subcluster)
        Lc_density_list = []
        for file in density_files:
            Lc_density_list.append(load2Col(file))
        return Lc_density_list

    def list_curve_shifts(self, subcluster):
        ''' 
        Returns a list of all shifts for this subcluster.
        '''
        coa = self.co_analysis
        indexes = self.list_curve_indexes(subcluster)
        shift = []
        i = 0
        while i < len(indexes):
            k = i+1
            while k < len(indexes):
                shift.append(round(coa.map_pair_to_shift(indexes[i], indexes[k]),4))
                k = k+1
            i = i+1
        return shift

    def list_initial_shifts(self, subcluster):
        '''
        Returns a list of the initial shifts to apply to align curves
        in the specified subcluster.
        '''
        shift_list = self.list_curve_shifts(subcluster)
        N = self.get_cluster_size(subcluster)
        initial_shifts = [0] + shift_list[0:N-1]
        return initial_shifts

    def plot_subcluster(self, subcluster, max_x = 150, alpha = 0.5, average=False):
        import matplotlib.pyplot as plt
        from numpy import average
        Lc_density_list = self.get_Lc_density_arrays(subcluster)
        initial_shifts = self.list_initial_shifts(subcluster)
        av_shift = average(initial_shifts)
        co = self.get_cluster_coincidence(subcluster)
        i = 0
        for Lc_density in Lc_density_list:
            plt.plot(Lc_density[0] + initial_shifts[i] - av_shift, Lc_density[1], 'k-', alpha=alpha)
            i = i+1

        plt.xlim(0,max_x)
        
        max_y = plt.axis()[3]
        plt.yticks([])
        plt.title("%s of %d curves ($\Gamma\geq$ %.3f)" % \
             (subcluster, self.get_cluster_size(subcluster), co))
        
        #plt.text(max_x/1.5, max_y/1.5, '$\Gamma\geq$ %.3f' % co, size = 10, color='blue')
        # plt.xlabel("Contour Length (nm)")
        # plt.ylabel("Density")
        # plt.title("Subcluster %d (%d %s) of %d curves" % \
        #     (self.get_cluster_node(subcluster).get_id(), self.get_cluster_node().get_id(), \
        #         subcluster, self.get_cluster_size(subcluster)))
        return plt.gcf()

    def plot_subcluster_noshift(self, subcluster, max_x = 150, alpha = 0.5):
        import matplotlib.pyplot as plt
        from numpy import average
        Lc_density_list = self.get_Lc_density_arrays(subcluster)
        i = 0
        for Lc_density in Lc_density_list:
            plt.plot(Lc_density[0], Lc_density[1], 'k-', alpha=alpha)
            i = i+1
        plt.xlim(0,max_x)
        plt.xlabel("Contour Length (nm)")
        plt.ylabel("Density")
        plt.title("Subcluster %d (%d %s) of %d curves" % \
            (self.get_cluster_node(subcluster).get_id(), self.get_cluster_node().get_id(), \
                subcluster, self.get_cluster_size(subcluster)))
        return plt.gcf()

    def plot_both_subclusters(self, max_x = 150, alpha=0.5):
        from matplotlib.pyplot import subplot, figsize
        #figsize(12,4)
        subplot(121)
        self.plot_subcluster('left', max_x = max_x, alpha=alpha, average=False)
        subplot(122)
        self.plot_subcluster('right', max_x = max_x, alpha=alpha, average=False)

''' 
Extraneous functions not in any object.
--------------------------------------
'''

def load_coincidence(fileName):
    #checking these
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
    ''' Takes a file with two columns and return each as an array. '''
    from numpy import asarray
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
    '''
    Import all the filenames of a certain type.
    '''
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
    '''Wrapper on fcluster.'''
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

def any_pre_order(root, node):
    '''
    Takes a root node of a whole tree and any specific node therein
    and returns a list of all the leaves in that particular node.
    '''
    import numpy as np
    n = root.count
    func=(lambda x: x.id)
    
    curNode = [None] * (2 * n)
    lvisited = np.zeros((2 * n,), dtype=bool)
    rvisited = np.zeros((2 * n,), dtype=bool)
    curNode[0] = node
    k = 0
    preorder = []
    while k >= 0:
        nd = curNode[k]
        ndid = nd.id
        if nd.is_leaf():
            preorder.append(func(nd))
            k = k - 1
        else:
            if not lvisited[ndid]:
                curNode[k + 1] = nd.left
                lvisited[ndid] = True
                k = k + 1
            elif not rvisited[ndid]:
                curNode[k + 1] = nd.right
                rvisited[ndid] = True
                k = k + 1
            # If we've visited the left and right of this non-leaf
            # node already, go up in the tree.
            else:
                k = k - 1
    return preorder

def get_cluster_node(root, left, right, curves):
    ''' 
    Takes a root node, the starting left and right nodes, and the curve array 
    you are trying to find.
    Returns the corresponding node index in the tree, and the two contained
    clusters (left and right) by curve index lists.
    '''
    from numpy import asarray, sort
    left_list = sort(asarray(any_pre_order(root, left)))
    right_list = sort(asarray(any_pre_order(root, right)))
    if len(left_list) == len(curves) and all(curves == left_list):
        return left
    elif len(right_list) == len(curves) and all(curves == right_list):
        return right
    else:
        if len(left_list) > len(curves) and any(curves[0] == left_list):
            return get_cluster_node(root, left.left, left.right, curves)
        elif len(right_list) > len(curves) and any(curves[0] == right_list):
            return get_cluster_node(root, right.left, right.right, curves)
        else:
            return "Error?"

''' --------------------- '''