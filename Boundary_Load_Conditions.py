#----------------------------------Displacement Driven------------------------------------------#
#----------- Following code is working for Coupling with connectivity matrix in progress-------------------#
import numpy as np 
import math
from Geometry import *


#Initializing global displcement vector
U_g_0=np.zeros(((nudof+nedof)*ncp,1))

#-------------Displacement Loading-------------#

#DispLoad1=[13,15,17]
#DispLoad2=[2,6]

#U_g_0[DispLoad1] = 0.1

#U_g_0[2][0] = 0.1
#U_g_0[5][0] = 0.1
#U_g_0[6][0] = 0.1
#U_g_0[7][0] = 0.1
#-------------Boundary Conditions--------------#
#---------------Displacement BCS---------------#
DispBcsFixed=[0,1,3,5,6,12]

U_g_0[DispBcsFixed] = 0

# U_g_0[0][0]     = 0
# U_g_0[1][0]     = 0
# U_g_0[3][0]     = 0
# U_g_0[4][0]     = 0
#----------------Electrical BCS----------------#
ElecBcsFixed=[18,21,24]

U_g_0[ElecBcsFixed] = 0

# U_g_0[8][0]     = 0
# U_g_0[10][0]    = 0

#----------------Electrical Loading----------------#

ElecLoad = [20,23,26]
U_g_0[ElecLoad] = 2e10

#U_g_0[9][0]     = 0
#U_g_0[11][0]    = 0


#-------------Replace this every time-----------------------#
BCS=np.sort(np.concatenate((DispBcsFixed,ElecBcsFixed,ElecLoad)))
#BCS = [0,1,3,4,5,7,8,10]
print(BCS)

print(U_g_0)
