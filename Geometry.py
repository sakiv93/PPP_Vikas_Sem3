#----------------------------------Displacement Driven------------------------------------------#
#------------------------------------------26 Feb-----------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from preprocessing_functions import *
from scipy import special
np.set_printoptions(threshold=np.inf)
#----------Manually change below variables--------------#

#------------------------------------Degree of the curve----------------------------------------#
p=2 #Degree of the curve in xi direction
q=2 #Degree of the curve in eta direction

#-----------------------------------Dimentions of 2D Plate--------------------------------------#
Thick   = 1.0 #Thickness of the plate
Length  = 10.0
Height  = 10.0

#------------------Number of Control points in xi and eta direction-----------------------------#
ncpxi    = 3 # No.of control points in xi  direction
ncpeta   = 3 # No.of control points in eta direction

class Switcher(object):
          def indirect(self,i):
                   method_name='Gauss_'+str(i)
                   method=getattr(self,method_name,lambda :'Invalid')
                   return method()
          def Gauss_0(self):
                   return np.array([[0],[2]])
          def Gauss_1(self):
                   return np.array([[0.5773,-0.5773],[1.0,1.0]])
          def Gauss_2(self):
                   return np.array([[0.7746,-0.7746,0.0],[0.5556,0.5556,0.8889]])
          def Gauss_3(self):
                   return np.array([[0.8611,-0.8611,0.3399,-0.3399],[0.3478, 0.3478,0.6521,0.6521]])

def GaussPoints(p,q):
    GaussMatrix=np.zeros(((p+1)*(q+1),3))
    s=Switcher()
    GaussPoints_xi  = s.indirect(p)
    GaussPoints_eta = s.indirect(q)
    print(GaussPoints_xi)
    print(GaussPoints_eta)
    k=0
    for i in range(p+1):
        for j in range(q+1):
            GaussMatrix[k,0]= GaussPoints_xi[0,i]
            GaussMatrix[k,1]= GaussPoints_eta[0,j]
            GaussMatrix[k,2]= GaussPoints_xi[1,i] * GaussPoints_eta[1,j]
            k+=1       
    return GaussMatrix

GPs_Ws = GaussPoints(p,q)

# if(p==1 and q==1):
#     #GPs_Ws = np.array([[-0.57735,-0.57735,1],[0.57735,-0.57735,1],[0.57735,0.57735,1],[-0.57735,0.57735,1]])
#     GPs_Ws = np.array([[-1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),1/(3**(1/2)),1],[-1/(3**(1/2)),1/(3**(1/2)),1]])

# else:
#     GPs_Ws = np.array([[-1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),-1/(3**(1/2)),1],[1/(3**(1/2)),1/(3**(1/2)),1],[-1/(3**(1/2)),1/(3**(1/2)),1]])

#     # GPs_Ws = np.array([[0.7746,0.7746,0.2743],[0.7746,-0.7746,0.2743],[0.7746,0.,0.4390],[-0.7746,0.7746,0.2743],
#     #                     [-0.7746,-0.7746,0.2743],[-0.7746,0.,0.4390],[0.,0.7746,0.4390],[0.,-0.7746,0.4390],[0.,0.,0.7023]])

p_ord=p+1 #order of the curve in xi direction
q_ord=q+1 #order of the curve in eta direction

d=1 #Derivatives need upto (k+l=d) Only first derivatives are needed for Analysis
nudof = 2 #Number of displacement degrees of freedom on each node
nedof = 1 #Number of electrical degrees of freedom on each node

#----------Need not alter below variables--------------#
#-----------------------Calculating Control point matrix----------------------------------#
P_W =np.zeros((ncpeta,ncpxi,4)) 
for j in range(ncpeta):
    for i in range(ncpxi):
        P_W[j,i,0] = (Length/(ncpxi-1))*i
        P_W[j,i,1] = (Height/(ncpeta-1))*j
        P_W[j,i,2] = 0
        P_W[j,i,3] = 1                     # Weights for respective control points
print(P_W)

# Input control point vector to element routine as a transpose
P_W_T = P_W.transpose((1,0,2))  #Here it is a 3D array (0,1,2) -- it is transposed to (1,0,2)

#******Check how to dicide points in x and y direction*****************
# ncpxi=np.shape(P_W)[1]  #No.of control points xi direction
# ncpeta=np.shape(P_W)[0] #No.of control points eta direction

U = KnotVector(p,ncpxi)
V = KnotVector(q,ncpeta)
print(U)
print(V)


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

#****** Algorithm from IGA Simplified paper ********#

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