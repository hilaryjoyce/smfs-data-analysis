from Coincidence.coincidence import *
#from Transform.Lc_transform import *
#from Transform.Plateau_transform import *
from Density.kernel_density import *

super_folder = '../Data/DNA_tiny/Transform_analysis/Parameters_SD0.017/'
#saveKDE(super_folder, covfac = 0.3, delta = 0.001, min_x = -1, max_x = 1)
saveCoincidence(super_folder, shift=0)