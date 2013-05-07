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

    def ordered_clusters(self):
        from numpy import asarray, argsort
        sizes = asarray(self.cluster_sizes())
        clusters = asarray(self.cluster_list())
        args = argsort(sizes)[::-1]
        return sizes[args], clusters[args]

    def curves_by_cluster(self):
        clusters = self.cluster_list()
        flat = self.flat_cluster()
        N = len(self.T)
        curve_list = [] * self.number_of_clusters()
        for i in range(N):
            curve_list[flat[i]+1] +=

class Cluster:

    def __init__(self, cluster_id, flat_cluster):
        
        self.id = cluster_id







