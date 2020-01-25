#-------------------------------Displacement Driven----------------------------------------#
#----------- Following code is working for Displacement Driven Coupling -------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
from Material_Properties import *

def materialRoutine(epsilon,electric_field, T_m):
    mu = youngs_modulus/(2*(1+poissons_ratio))
    lamda = (poissons_ratio*youngs_modulus)/((1+poissons_ratio)*(1-2*poissons_ratio))
    sigma_y = yield_stress
    # Calculating Material tangent stiffness matrix
    Ct = np.array([[2*mu+lamda,lamda,0],[lamda,2*mu+lamda,0],[0,0,mu]])
    #print('Ct',Ct)
    #Ct = np.array([[139000,74280,0],[74280,115400,0],[0,0,25640]])

    sigma  = np.matmul(Ct,epsilon)-np.matmul(np.transpose(e),electric_field)
    Electrical_Displacement = np.matmul(e,epsilon)+np.matmul(k,electric_field)


    # returning Material tangent stiffness matrix, Updated stress to element routine,
    # Piezoelectric constant matrix, Dielectric material constants
    return Ct, e, k, sigma, Electrical_Displacement



# C=materialRoutine([[0],[0],[0]], 1)
# print(C)