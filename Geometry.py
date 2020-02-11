#----------------------------------Displacement Driven------------------------------------------#
#---------------------------Code Works for any degree of the NURBS Curve------------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from preprocessing_functions import *
from scipy import special
np.set_printoptions(threshold=np.inf)
#----------Manually change below variables--------------#

Thick   = 1.0 #Thickness of the plate
Length  = 10.0
Height  = 10.0

#------------------Number of elements in xi and eta direction-----------------------------#
nexi    = 1 # No.of elements in xi  direction
neeta   = 1 # No.of elements in eta direction

# GPs_Ws = np.array([[-0.57735,-0.57735,1],[0.57735,-0.57735,1],[0.57735,0.57735,1],[-0.57735,0.57735,1]])
GPs_Ws = np.array([[0.7746,0.7746,0.2743],[0.7746,-0.7746,0.2743],[0.7746,0.,0.4390],[-0.7746,0.7746,0.2743],
                    [-0.7746,-0.7746,0.2743],[-0.7746,0.,0.4390],[0.,0.7746,0.4390],[0.,-0.7746,0.4390],[0.,0.,0.7023]])


p=3 #Degree of the curve in xi direction
q=3 #Degree of the curve in eta direction

p_ord=p+1 #order of the curve in xi direction
q_ord=q+1 #order of the curve in eta direction

d=1 #Derivatives need upto (k+l=d) Only first derivatives are needed for Analysis
nudof = 2 #Number of displacement degrees of freedom on each node
nedof = 1 #Number of electrical degrees of freedom on each node

#----------Need not alter below variables--------------#
#-----------------------Calculating Control point matrix----------------------------------#
# P_W=np.array([[[0.,0.,0.,1.],[5.,0.,0.,1.],[10.,0.,0.,1.]],
#              [[0.,5.,0.,1.],[5.,5.,0.,1.],[10.,5.,0.,1.]],
#              [[0.,10.,0.,1.],[5.,10.,0.,1.],[10.,10.,0.,1.]]])

# P_W=np.array([[[0,0,0,1],[10,0,0,1]],
#              [[0,5,0,1],[10,5,0,1]],
#              [[0,10,0,1],[10,10,0,1]]])

P_W=np.array([[[0.,0.,0.,1.],[2.5,0.,0.,1.],[5.,0.,0.,1.],[10.,0.,0.,1.]],
             [[0.,2.5,0.,1.],[2.5,2.5,0.,1.],[5.,2.5,0.,1.],[10.,2.5,0.,1.]],
             [[0.,5.,0.,1.],[2.5,5.,0.,1.],[5.,5.,0.,1.],[10.,5.,0.,1.]],
             [[0.,10.,0.,1.],[2.5,10.,0.,1.],[5.,10.,0.,1.],[10.,10.,0.,1.]]])

# Input control point vector to element routine as a transpose
P_W_T = P_W.transpose((1,0,2))  #Here it is a 3D array (0,1,2) -- it is transposed to (1,0,2)

#******Check how to dicide points in x and y direction*****************
ncpxi=np.shape(P_W)[1]  #No.of control points xi direction
ncpeta=np.shape(P_W)[0] #No.of control points eta direction

U = KnotVector(p,ncpxi)
V = KnotVector(q,ncpeta)
#U = np.array([0., 0., 0., 1, 2., 2., 2.])  #Knot vector in x direction
#V = np.array([0., 0., 0., 1, 2., 2., 2.])   #Knot vector in y direction
print('KnotVectorU',U)
print('KnotVectorV',V)


n=(np.size(U)-1)-p-1
m=(np.size(V)-1)-q-1
P=P_W_T[:,:,0:3]
W=P_W_T[:,:,3]

#**********************************************************************
ncp = ncpxi*ncpeta #Total number of control points

nel= (ncpxi-p)*(ncpeta-q)
necp= (p+1)*(q+1)   #Total number of control points per element
necpxi  = (p+1)     #Number of control points per element along xi  direction
necpeta = (q+1)     #Number of control points per element along eta direction

#------------------------Knot connectivity-----------------------------#

uniqueU = np.unique(U)
uniqueV = np.unique(V)

nelU = ncpxi-p
nelV = ncpeta-q
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
knotConnectivity = knotConnectivity -1
knotConnectivity = knotConnectivity.astype(int)
print('knotConnectivity',knotConnectivity)