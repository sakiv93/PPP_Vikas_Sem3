#-------------------------------Displacement Driven----------------------------------------#
#----------- Following code is working for Displacement Driven Coupling -------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
#----------Manually change below variables--------------#

#GPs_Ws = np.array([[-0.57735,-0.57735,1],[0.57735,-0.57735,1],[0.57735,0.57735,1],[-0.57735,0.57735,1]])
GPs_Ws = np.array([[-1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),1/(3**(1/2)),1],[-1/(3**(1/2)),1/(3**(1/2)),1]])

p=1 #Degree of the curve in x direction
q=1 #Degree of the curve in y direction

U = np.array([0., 0., 1.,1.])  #Knot vector in x direction
V = np.array([0., 0., 1.,1.])  #Knot vector in y direction

P_W=np.array([[[0,0,0,1],[10,0,0,1]],
          [[0,10,0,1],[10,10,0,1]]])

P_W_T = P_W.transpose((1,0,2))  #Here it is a 3D array (0,1,2) -- it is transposed to (1,0,2)
# Input control point vector to element routine as a transpose #
# P_W=np.array([[[0,0,0,1],[0,1,0,1]],
#             [[1,0,0,1],[1,1,0,1]]])




d=1 #Derivatives need upto (k+l=d)
nudof = 2 #Number of displacement degrees of freedom on each node
nedof = 1 #Number of electrical degrees of freedom on each node


#----------Need not alter below variables--------------#
n=(np.size(U)-1)-p-1
m=(np.size(V)-1)-q-1
P=P_W_T[:,:,0:3]
W=P_W_T[:,:,3]
#******Check how to dicide points in x and y direction*****************
ncpxi=np.shape(P_W)[1]  #No.of control points xi direction
ncpeta=np.shape(P_W)[0] #No.of control points eta direction
#**********************************************************************
#nel=(ncpxi-p-1)*(ncpeta-q-1)
ncp = ncpxi*ncpeta #Total number of control points
#necp = (p+2)*(q+2) #Total number of control points per element
nel=1
necp=4 #Total number of control points per element

Thick=1.0 #Thickness of the plate