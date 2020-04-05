#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#
import numpy as np 
import math
from Geometry import *

#Initializing global displcement vector
U_g_0=np.zeros(((nudof+nedof)*ncp,1))

# The following code isolates Degrees of freedom numbers of nodes 
# for easily handling bounday conditions (BCS) definition. 

#-------------Bottom nodes Y displacement Dof----------------#
Bottom_nodes_u = np.zeros(ncpxi)
for BNu in range(ncpxi):
    Bottom_nodes_u[BNu] =2*BNu+1
Bottom_nodes_u= Bottom_nodes_u.astype(int)

#-------------Left nodes X displacement Dof------------------#
Left_nodes_u = np.zeros(ncpeta)
for LNu in range(ncpeta):
    Left_nodes_u[LNu] = (ncpxi*2)*LNu
Left_nodes_u= Left_nodes_u.astype(int)

#-------------Right nodes X displacement  Dof------------------#
Right_nodes_u = np.zeros(ncpeta)
for RNu in range(ncpeta):
    Right_nodes_u[RNu] = ((ncpxi*2)*RNu) + (ncpxi-1)*2
Right_nodes_u= Right_nodes_u.astype(int)

#-------------Top nodes Y displacement Dof----------------#
Top_nodes_u = np.zeros(ncpxi)
for TNu in range(ncpxi):
    Top_nodes_u[TNu] =((ncpxi*2) * (ncpeta-1) + 1) + (TNu*2)
Top_nodes_u= Top_nodes_u.astype(int)

#-------------Left nodes Electric Dof------------------#
Left_nodes_e = np.zeros(ncpeta)
for LNe in range(ncpeta):
    Left_nodes_e[LNe] = ncpxi*ncpeta*2 + ncpxi*LNe
Left_nodes_e= Left_nodes_e.astype(int)

#-------------Right nodes Electric Dof------------------#
Right_nodes_e = np.zeros(ncpeta)
for RNe in range(ncpeta):
    Right_nodes_e[RNe] = (ncpxi*ncpeta*2 -1) + ncpxi*(RNe+1)
Right_nodes_e= Right_nodes_e.astype(int)

#---------------Top nodes Electric Dof-------------------#
Top_nodes_e = np.zeros(ncpxi)
for TNe in range(ncpxi):
    Top_nodes_e[TNe] = ncpxi*ncpeta*2 + ncpxi*(ncpeta-1) +TNe
Top_nodes_e= Top_nodes_e.astype(int)

#-------------Bottom nodes Electric Dof------------------#
Bottom_nodes_e = np.zeros(ncpxi)
for BNe in range(ncpxi):
    Bottom_nodes_e[BNe] = ncpxi*ncpeta*2 + BNe
Bottom_nodes_e= Bottom_nodes_e.astype(int)

#---------------------Boundary Conditions definition-------------------------#

#-------------Displacement Loading-------------#
#U_g_0[Right_nodes_u]  = 1e-4
#U_g_0[Top_nodes_u]    = 0.2

#---------------Fixed Displacement BCS---------------#
U_g_0[Bottom_nodes_u] = 0
U_g_0[Left_nodes_u]   = 0

#----------------Grounded Electrical BCS----------------#
U_g_0[Top_nodes_e]    = 0

#----------------Electrical Loading----------------#
U_g_0[Bottom_nodes_e]  = 100

#-------------Replace this every time for defining BCS in main program-----------------------#
# Or direct node numbers can be given where boundary conditions for DOF are given 
# Format [2,10,3] (Example : DOF numbers inside an list)
BCS=np.sort(np.concatenate((Bottom_nodes_u,Left_nodes_u,Top_nodes_e,Bottom_nodes_e)))
print(BCS)
