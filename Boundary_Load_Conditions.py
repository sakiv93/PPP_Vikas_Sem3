#----------------------------------Displacement Driven------------------------------------------#
#----------------------------One Program for any degree curve-----------------------------------#
import numpy as np 
import math
from Geometry import *

#Initializing global displcement vector
U_g_0=np.zeros(((nudof+nedof)*ncp,1))

#-------------Bottom nodes displacement y Dof----------------#
Bottom_nodes_u = np.zeros(ncpxi)
for BNu in range(ncpxi):
    Bottom_nodes_u[BNu] =2*BNu+1
Bottom_nodes_u= Bottom_nodes_u.astype(int)
print('Bottom_nodes_u',Bottom_nodes_u)

#-------------Left nodes displacement x Dof------------------#
Left_nodes_u = np.zeros(ncpeta)
for LNu in range(ncpeta):
    Left_nodes_u[LNu] = (ncpxi*2)*LNu
Left_nodes_u= Left_nodes_u.astype(int)
print(Left_nodes_u)

#-------------Right nodes displacement x Dof------------------#
Right_nodes_u = np.zeros(ncpeta)
for RNu in range(ncpeta):
    Right_nodes_u[RNu] = ((ncpxi*2)*RNu) + (ncpxi-1)*2
Right_nodes_u= Right_nodes_u.astype(int)
print(Right_nodes_u)

#-------------Top nodes displacement y Dof----------------#
Top_nodes_u = np.zeros(ncpxi)
for TNu in range(ncpxi):
    Top_nodes_u[TNu] =((ncpxi*2) * (ncpeta-1) + 1) + (TNu*2)
Top_nodes_u= Top_nodes_u.astype(int)
print(Top_nodes_u)

#-------------Left nodes Electric Dof------------------#
Left_nodes_e = np.zeros(ncpeta)
for LNe in range(ncpeta):
    Left_nodes_e[LNe] = ncpxi*ncpeta*2 + ncpxi*LNe
Left_nodes_e= Left_nodes_e.astype(int)
print(Left_nodes_e)

#-------------Right nodes Electric Dof------------------#
Right_nodes_e = np.zeros(ncpeta)
for RNe in range(ncpeta):
    Right_nodes_e[RNe] = (ncpxi*ncpeta*2 -1) + ncpxi*(RNe+1)
Right_nodes_e= Right_nodes_e.astype(int)
print(Right_nodes_e)

#---------------------Boundary Conditions-------------------------#

#-------------Displacement Loading-------------#
U_g_0[Right_nodes_u]  = 0.1
U_g_0[Top_nodes_u]    = 0.2

#---------------Displacement BCS---------------#
U_g_0[Bottom_nodes_u] = 0
U_g_0[Left_nodes_u]   = 0

#----------------Electrical BCS----------------#
U_g_0[Left_nodes_e]   = 0

#----------------Electrical Loading----------------#
#U_g_0[Right_nodes_e]  = 500000000

#-------------Replace this every time-----------------------#
BCS=np.sort(np.concatenate((Right_nodes_u,Top_nodes_u,Bottom_nodes_u,Left_nodes_u,Left_nodes_e)))
#print(BCS)
#print(U_g_0)
