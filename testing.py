from Coincidence.coincidence import *
from Transform.Lc_transform import *
from Transform.Plateau_transform import *
from Density.kernel_density import *

super_folder = '../Data/DNA_tiny/Transform_analysis/Parameters_SD0.017/'

saveKDE(super_folder, covfac = 0.3, delta = 0.001, min_x = -5.0, max_x = 5.0)