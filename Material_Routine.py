#----------------------------------Displacement Driven------------------------------------------#
#------------------------------------------26 Feb-----------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
from Material_Properties import *

def materialRoutine(epsilon,electric_field, T_m):
    """
    Input: 
        Input strain data(epsilon), Electric field data(electric_field) from element routine
    Process: 
        Material routine calculates Element tangent stiffness matrix 
        (Elastic (Ct), Piezoelectric (e) and dielectric constants (k))
        Stress (Sigma) and Electrical displacements (Electrical_Displacement) 
        are also calculated.
    Return: 
        The function returns Elastic, Piezoelectric, dielectric constant matrices, 
        Stress and Electrical displacements 
    """
    mu = youngs_modulus/(2*(1+poissons_ratio))
    lamda = (poissons_ratio*youngs_modulus)/((1+poissons_ratio)*(1-2*poissons_ratio))
    sigma_y = yield_stress

    # Calculating Material tangent stiffness matrix
    Ct = np.array([[2*mu+lamda,lamda,0],[lamda,2*mu+lamda,0],[0,0,mu]])

    # Calculating stress and electrical displacements
    sigma  = np.matmul(Ct,epsilon)-np.matmul(np.transpose(e),electric_field)
    Electrical_Displacement = np.matmul(e,epsilon)+np.matmul(k,electric_field)


    # Returning Elastic,Piezoelectric, Dielectric constants, Stress and Electrical Displacements
    return Ct, e, k, sigma, Electrical_Displacement