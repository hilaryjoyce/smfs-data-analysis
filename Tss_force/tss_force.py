''' tss_force.py contains the output functions to save tss-force data. 
    Parameters:
        - folder of raw data
        - type = ['RDF' or 'ADR']
        - spring constant for ADR type

    Currently only written for 'RDF' type.
'''

def tss_force_list(folder, type, k = 0):
    ''' Takes a folder of raw data and a type ('RDF' or 'ADR') and
        returns a list of paired tss-force arrays and a second list of the
        curve numbers for each file.'''
    from glob import glob
    if type == 'RDF':
        import RDF_pulling as pull
        
        RDF_folder = '%sRDF_text/' % folder
        RDF_files = glob('%s*.txt' % RDF_folder)

        tf_list = []
        curve_number_list = []
        for file in RDF_files:
            curve_number_list.append(pull.curveNumFinder(file))
            tf_list.append(pull.tss_force(file))

    if type == 'ADR':
        import ADR_pulling as pull
        
        ADR_folder = '%sADR_text/' % folder
        ADR_files = glob('%s*.txt' % ADR_folder)
        tf_list = []
        curve_number_list = []

        for file in ADR_files:
            curve_number_list.append(pull.curveNumFinder(file))
            tf_list.append(pull.tss_force(file, k=k))

    return tf_list, curve_number_list

def saveFile2(path, col1, col2, col1_name, col2_name):
    '''Saves a file with 2 columns of data.'''
    file = open(path, 'w')
    file.write("%s\t%s\n" % (col1_name, col2_name))
    k=0
    while k<len(col1):
        el = '%3g\t%3g\n' % (col1[k], col2[k])
        file.write(el)
        k=k+1
    file.close

def saveTssForce(folder, type, k = 0):
    ''' Save a .txt file of tss and force to Tss_force_data folder.
        Should only be executed once for all files in
        a given experiment or it's rather redundant.'''
    import os
    if not os.path.isdir("%sTss_Force_text/" % folder):
        os.mkdir("%sTss_Force_text/" % folder)

    tf_list, curve_number_list = tss_force_list(folder, type, k=k)
    i = 0
    while i < len(curve_number_list):
        filename = "%sTss_Force_text/tss_force_%s.txt" % (folder, curve_number_list[i])
        saveFile2(filename, tf_list[i][0], tf_list[i][1], "Tss (nm)", "Force (pN)")
        i = i+1

