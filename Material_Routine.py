#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#

#-----------------------------------------------------------------------------------------------#
#                               List of variables used                                          #
#           Ct- Elastic constants                                                               #
#           e - Piezoelectric constants                                                         #
#           k - dielectric constants                                                            #
#           sigma - stress                                                                      #
#           Electrical_Displacement - Electrical Displacement                                   #
#           epsilon - strain values                                                             #
#           electric_field - electric field values                                              #
#-----------------------------------------------------------------------------------------------#

import numpy as np
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

    # Calculating stress and electrical displacements
    sigma  = np.matmul(Ct,epsilon)-np.matmul(np.transpose(e),electric_field)
    print('np.matmul(np.transpose(e),electric_field)',np.matmul(np.transpose(e),electric_field))
    print('np.matmul(Ct,epsilon)',np.matmul(Ct,epsilon))
    Electrical_Displacement = np.matmul(e,epsilon)+np.matmul(k,electric_field)


    # Returning Elastic,Piezoelectric, Dielectric constants, Stress and Electrical Displacements
    return Ct, e, k, sigma, Electrical_Displacement