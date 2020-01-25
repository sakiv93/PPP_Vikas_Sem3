#-------------------------------Displacement Driven----------------------------------------#
#----------- Following code is working for Displacement Driven Coupling -------------------#
import numpy as np 
import math
from Geometry import *


#Initializing global displcement vector
U_g_0=np.zeros(((nudof+nedof)*necp,1))

#-------------Displacement Loading-------------#
#U_g_0[5][0] = 0.1
#U_g_0[7][0] = 0.1
U_g_0[5][0] = 0.1
U_g_0[7][0] = 0.1
#-------------Boundary Conditions--------------#
#---------------Displacement BCS---------------#
U_g_0[0][0]     = 0
U_g_0[1][0]     = 0
U_g_0[3][0]     = 0
U_g_0[4][0]     = 0
#----------------Electrical BCS----------------#
U_g_0[8][0]     = 0
U_g_0[10][0]    = 0

#----------------Electrical Loading----------------#
#U_g_0[9][0]     = 1000
#U_g_0[11][0]    = 1000


BCS=[0,1,3,4,5,7,8,10]
