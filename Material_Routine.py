import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
from Material_Properties import *

def materialRoutine(epsilon, T_m):
    mu = youngs_modulus/(2*(1+poissons_ratio))
    lamda = (poissons_ratio*youngs_modulus)/((1+poissons_ratio)*(1-2*poissons_ratio))
    sigma_y = yield_stress
    # Calculating Material tangent stiffness matrix
    Ct = np.array([[2*mu+lamda,lamda,0],[lamda,2*mu+lamda,0],[0,0,mu]]) 
    sigma  = np.matmul(Ct,epsilon) 
    # returning Material tangent stiffness matrix and Updated stress to element routine
    return Ct, sigma

# C=materialRoutine([[0],[0],[0]], 1)
# print(C)