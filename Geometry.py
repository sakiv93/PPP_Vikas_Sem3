#----------------------------------Displacement Driven------------------------------------------#
#----------- Following code is working for Coupling with connectivity matrix in progress-------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
#----------Manually change below variables--------------#

#GPs_Ws = np.array([[-0.57735,-0.57735,1],[0.57735,-0.57735,1],[0.57735,0.57735,1],[-0.57735,0.57735,1]])
GPs_Ws = np.array([[-1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),1/(3**(1/2)),1],[-1/(3**(1/2)),1/(3**(1/2)),1]])

p=1 #Degree of the curve in xi direction
q=1 #Degree of the curve in eta direction

p_ord=p+1 #order of the curve in xi direction
q_ord=q+1 #order of the curve in eta direction

U = np.array([0.,0.,1.,2.,2.])  #Knot vector in x direction
V = np.array([0.,0.,1.,2.,2.])  #Knot vector in y direction

P_W=np.array([[[0,0,0,1],[5,0,0,1],[10,0,0,1]],
          [[0,5,0,1],[5,5,0,1],[10,5,0,1]],
          [[0,10,0,1],[5,10,0,1],[10,10,0,1]]])

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
#### nel=1
nel= (ncpxi-p_ord)*(ncpeta-q_ord)
necp= (p_ord+1)*(q_ord+1) #Total number of control points per element

Thick=1.0 #Thickness of the plate


#------------------------Knot connectivity-----------------------------#

uniqueU = np.unique(U)
uniqueV = np.unique(V)
#print(np.unique(U))
#print(np.unique(V))
nelU = ncpxi-p_ord
nelV = ncpeta-q_ord
nel = nelU*nelV
knotConnectivity = np.zeros((nel,2)) 
Span_U = np.zeros((nelU,2)) # Span_U have rows equal to number of elements in xi direction
Span_V = np.zeros((nelV,2)) # Span_V have rows equal to number of elements in eta direction
count=0
for i in range(0,nelV):
    for j in range(0,nelU):
        knotConnectivity[count,:] = [j+1,i+1] 
        # First and second coloumns of 'knotConnectivity' store the global number of knot span ranges
        # or row number of Span_U and Span_V arrays
        Span_U[j,:]= [uniqueU[j],uniqueU[j+1]]
        Span_V[i,:]= [uniqueV[i],uniqueV[i+1]]
        count=count+1