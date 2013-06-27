'''
split_all_coincidence.py
Takes an all_coincidence_report.txt IN a Coincidence_max_x folder already
and returns a folder + coincidence report for each shift.
'''

def split_all_coincidence(param_folder, type='density', max_x=300):
    '''
    Takes an all_coincidence_report.txt IN a Coincidence_max_x folder already
    and returns a folder + coincidence report for each shift.
    '''

    from glob import glob
    import os 
    
    if type == 'density':
        co_folder = '%sCoincidence_max%d/' % (param_folder, max_x)
    else:
        co_folder = '%sCoincidence/' % param_folder
    co_file = glob('%sall*.txt' % co_folder)[0]

    file = open(co_file, 'r')
    header = file.readline().split()

    # deal with the header
    sg_header = header[2:]
    shift_list = []
    i = 0
    while i < len(sg_header):
        try:
            shift_list.append(int(sg_header[i][1:]))
        except:
            shift_list.append(sg_header[i][1:])
        i = i+2
    print shift_list
    N = len(shift_list)

    shift_data = []
    i = 0
    while i < N:
        shift_data.append([[],[]])
        i = i+1
        
    # break up the rows
    data = []
    curve1_list = []
    curve2_list = []
    for line in file:
        row = line.split()
        curve1_list.append(row[0])
        curve2_list.append(row[1])
        i = 0 
        while i < N:
            shift_data[i][0].append(float(row[i*2+2]))
            shift_data[i][1].append(float(row[i*2+3]))
            i = i+1
    file.close()

    shift_folders = []
    for shift in shift_list:
        if shift == 'No':
            savefolder = "NoShift/"
        else:
            savefolder = "Shift_%d/" % int(shift)
        
        if not os.path.isdir(co_folder+savefolder):
            os.mkdir(co_folder+savefolder)
        
        shift_folders.append(co_folder+savefolder)

    i = 0
    for data in shift_data:
        with open('%scoincidence_report.txt' % shift_folders[i], 'w') as file:
            file.write("# c1\tc2\tGamma\tBest shift (nm)\n")
            k = 0
            while k < len(data[0]):
                line = '%s\t%s\t%.3f\t%5.2f\n' % (curve1_list[k], curve2_list[k], data[0][k], data[1][k])
                file.write(line)
                k = k+1
        i = i+1

    return None

def split_one_coincidence(param_folder, id=0, type='density', max_x=300):
    '''
    Takes an all_coincidence_report.txt IN a Coincidence_max_x folder already
    and returns a folder + coincidence report for each shift.
    '''

    from glob import glob
    import os 
    import sys
    
    if type == 'density':
        co_folder = '%sCoincidence_max%d/' % (param_folder, max_x)
    else:
        co_folder = '%sCoincidence/' % param_folder
    co_file = glob('%sall*.txt' % co_folder)[0]

    file = open(co_file, 'r')
    header = file.readline().split()

    # deal with the header
    sg_header = header[2:]
    shift_list = []
    i = 0
    while i < len(sg_header):
        try:
            shift_list.append(int(sg_header[i][1:]))
        except:
            shift_list.append(sg_header[i][1:])
        i = i+2
    print shift_list
    N = len(shift_list)

    if id >= N:
        print "invalid id"
        sys.exit()

    i = 0
    print "Splitting shift = %s" % str(shift_list[id])

    shift_data = []
        
    # break up the rows
    curve1_list = []
    curve2_list = []
    i = 0
    for line in file:
        row = line.split()
        curve1_list.append(row[0])
        curve2_list.append(row[1])
        shift_data.append([float(row[id*2+2]), float(row[id*2+3])])
        if i%100000 == 0:
            print i
        i = i+1
    file.close()

    shift = shift_list[id]
    if shift == 'No':
        savefolder = "NoShift/"
    else:
        savefolder = "Shift_%d/" % int(shift)
    
    if not os.path.isdir(co_folder+savefolder):
        os.mkdir(co_folder+savefolder)
    
    print co_folder+savefolder
    with open('%scoincidence_report.txt' % (co_folder+savefolder), 'w') as file:
        file.write("# c1\tc2\tGamma\tBest shift (nm)\n")
        k = 0
        while k < len(shift_data):
            line = '%s\t%s\t%.3f\t%5.2f\n' % (curve1_list[k], curve2_list[k], shift_data[k][0], shift_data[k][1])
            file.write(line)
            k = k+1
    return None
