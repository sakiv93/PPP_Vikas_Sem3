#-------------------------------Displacement Driven----------------------------------------#
#----------- Following code is working for Displacement Driven Coupling -------------------#
import numpy as np 
import math
from Geometry import *


#Initializing global displcement vector
U_g_0=np.zeros(((nudof+nedof)*necp,1))

#-------------Displacement Loading-------------#
U_g_0[2][0] = 8e-04
U_g_0[6][0] = 8e-04

#-------------Boundary Conditions--------------#
#---------------Displacement BCS---------------#
U_g_0[0][0]     = 0
U_g_0[1][0]     = 0
U_g_0[4][0]     = 0
U_g_0[5][0]     = 0
#----------------Electrical BCS----------------#
U_g_0[8][0]     = 0
U_g_0[10][0]    = 0


BCS=[0,1,2,4,5,6,8,10]
