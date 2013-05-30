# move sample files

import shutil

curved_files = "curved_sample_files.txt"
flat_files = "flat_sample_files.txt"

curved_file = open(curved_files, 'r')
curved_list = []
for line in curved_file:
    curved_list.append(line.strip('\n'))
curved_file.close()

flat_file = open(flat_files, 'r')
flat_list = []
for line in flat_file:
    flat_list.append(line.strip('\n'))
flat_file.close()

dest_curved = "../Data/aLa_curved_sample/RDF_text/"
dest_flat = "../Data/aLa_flat_sample/RDF_text/"
dest_all = "../Data/aLa_all_sample/RDF_text/"

for file in curved_list:
    shutil.copy(file, dest_all)
    shutil.copy(file, dest_curved)

for file in flat_list:
    shutil.copy(file, dest_all)
    shutil.copy(file, dest_flat)