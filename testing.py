from Tss_force_code.tss_force import *

folder = '../Data/DNA_tiny/'

tf_list, cnum_list = tss_force_list(folder, type='ADR', k=0.367)
saveTssForce(folder, type='ADR', k = 0.367)
print cnum_list
print tf_list[0][0]