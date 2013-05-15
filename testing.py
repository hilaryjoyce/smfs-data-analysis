#from Coincidence.coincidence import *
#from Transform.Lc_transform import *
#from Transform.Plateau_transform import *
#from Density.kernel_density import *
from Tss_force.tss_force import saveTssForce

super_folder = '../Data/aLa_flat_all/'
saveTssForce(super_folder, type='RDF')
#saveKDE(super_folder, covfac = 0.3, delta = 0.001, min_x = -1, max_x = 1)
