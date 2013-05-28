#run ala_all 

#from Transform.Lc_transform import saveLcList
#from Density.kernel_density import saveKDE
#from Tss_force.tss_force import saveTssForce
from Coincidence.coincidence import saveMultipleCoincidence

super_folder = "../Data/aLa_all/"
param_folder = '../Data/aLa_all/Transform_analysis/Parameters_p0.4_minTss5_minForce45/'

saveMultipleCoincidence(param_folder)