#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#                               List of variables used                                          #
#           p,q             - Degree of the curve in xi and eta direction respectively          #
#           ncpxi, ncpeta   - No.of control points in xi and eta direction respectively         #
#           Pw              - Control points matrix                                             #
#           XI, ETA         - Knot Vectors in xi and eta direction respectively                 #
#           ncp             - Total number of control points                                    #
#           nel             - Total number of elements                                          #
#-----------------------------------------------------------------------------------------------#
import numpy as np
from preprocessing_functions import *
from Input import *

#----------Manually change below variables--------------#

# #------------------------------------Degree of the curve----------------------------------------#
# p=2 #Degree of the curve in xi direction
# q=2 #Degree of the curve in eta direction

# #-----------------------------------Dimentions of 2D Plate--------------------------------------#
# Thick   = 1.0  #Thickness of the plate in mm
# Length  = 10.0 #Length of the plate    in mm
# Height  = 10.0 #Height of the plate    in mm

# #------------------Number of Control points in xi and eta direction-----------------------------#
# ncpxi    = 2 # No.of control points in xi  direction
# ncpeta   = 2 # No.of control points in eta direction

#----------Need not alter below variables--------------#

# A switch class is implemented inorder to select the gauss points and weights according with 
# degree of the curve in each direction.
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

GPs_Ws = GaussPoints(p,q) # Gauss points along with their respective weights

p_ord=p+1 #order of the curve in xi direction
q_ord=q+1 #order of the curve in eta direction

d=1         #Derivatives need upto (k+l=d) Only first derivatives are needed for Analysis
nudof = 2   #Number of displacement degrees of freedom on each node
nedof = 1   #Number of electrical degrees of freedom on each node

#-----------------------Calculating Control point matrix----------------------------------#
Pw =np.zeros((ncpeta,ncpxi,4)) 
for j in range(ncpeta):
    for i in range(ncpxi):
        Pw[j,i,0] = (Length/(ncpxi-1))*i
        Pw[j,i,1] = (Height/(ncpeta-1))*j
        Pw[j,i,2] = 0
        Pw[j,i,3] = 1                     # Weights for respective control points

print(Pw)

# Input control point vector to element routine as a transpose
Pw_T = Pw.transpose((1,0,2))  #Here it is a 3D array (0,1,2) -- it is transposed to (1,0,2)

XI = KnotVector(p,ncpxi)    # Knot Vector in xi  direction
ETA = KnotVector(q,ncpeta)  # Knot Vector in eta direction
print(XI)
print(ETA)

n=(np.size(XI)-1)-p-1
m=(np.size(ETA)-1)-q-1
P=Pw_T[:,:,0:3]
W=Pw_T[:,:,3]

#**********************************************************************
ncp = ncpxi*ncpeta #Total number of control points

nel= (ncpxi-p)*(ncpeta-q)
necp= (p+1)*(q+1)   #Total number of control points per element
necpxi  = (p+1)     #Number of control points per element along xi  direction
necpeta = (q+1)     #Number of control points per element along eta direction

#------------------------Knot connectivity-----------------------------#

uniqueXI  = np.unique(XI)
uniqueETA = np.unique(ETA)

nelXI  = ncpxi-p
nelETA = ncpeta-q
knotConnectivity = np.zeros((nel,2)) 
Span_XI  = np.zeros((nelXI,2))  # Span_XI have rows equal to number of elements in xi direction
Span_ETA = np.zeros((nelETA,2)) # Span_ETA have rows equal to number of elements in eta direction
count=0

#-----------------------------------Algorithm adopted from -------------------------------------# 
# Vishal Agrawal and Sachin S Gautam. Iga: A simplied introduction and implementation
# details for nite element users. Journal of The Institution of Engineers (India): Series C,
# 100(3):561{585, 2019.

for i in range(0,nelETA):
    for j in range(0,nelXI):
        knotConnectivity[count,:] = [j+1,i+1] 
        # First and second coloumns of 'knotConnectivity' store the global number of knot span ranges
        # or row number of Span_XI and Span_ETA arrays
        Span_XI[j,:]= [uniqueXI[j],uniqueXI[j+1]]
        Span_ETA[i,:]= [uniqueETA[i],uniqueETA[i+1]]
        count=count+1
knotConnectivity = knotConnectivity -1
knotConnectivity = knotConnectivity.astype(int)
print('knotConnectivity',knotConnectivity)
print('Span_XI',Span_XI)
print('Span_ETA',Span_ETA)
#------------------------------------------------------------------------------------------------# 