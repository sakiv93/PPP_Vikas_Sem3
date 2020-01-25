#-------------------------------Displacement Driven----------------------------------------#
#----------- Following code is working for Displacement Driven Coupling -------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
from Geometry import *
from Material_Routine import *
from preprocessing_functions import *

def Jacobian12(j,xi,eta,elU,elV):

    #------------Calculating J1-------------------------#

    dxi_dximas = 0.5*(elU[1]-elU[0])
    deta_detamas = 0.5*(elV[1]-elV[0]) 
    J2det = dxi_dximas*deta_detamas

    #------------Calculating J2-------------------------#

    #--------Evaluation of NURBS Surface basis functions derivatives---------#
    Aders,wders = SurfaceDerivsAlgAuv(n,p,U,m,q,V,P,W,xi,eta,d)
    #print('Aders:',Aders)
    #print('wders:',wders)
    dRx = RatSurfaceDerivs(Aders,wders,d)
    dRx_dxi  =  dRx[1][0][0]
    dRx_deta =  dRx[0][1][0]
    dRy_dxi  =  dRx[1][0][1]
    dRy_deta =  dRx[0][1][1]
    J1 = np.array([[dRx_dxi,dRx_deta],
                [dRy_dxi,dRy_deta]])
    #print('J1 Matrix : ',J1)
    J1det = (dRx_dxi*dRy_deta)-(dRx_deta*dRy_dxi)
    #print(J1det)
    #print('Determinat of J1 : ',J1det)
    J1inv = np.linalg.inv(J1)
    #print('J1inv : ',J1inv)
    return J1det,J2det

def B_matrix(xi,eta,Bu,Be):

    #---------------NURBS Basis Functions Derivatives wrt x and y-------------#
    uspan = FindSpan(n,p,xi,U)
    vspan = FindSpan(m,q,eta,V)
    # Indices of non vanishing basis functions
    NVu = np.arange(uspan-p,uspan+1,1) #Thoroughly checked. Can crosscheck again
    NVv = np.arange(vspan-p,vspan+1,1)
    Aders,wders = SurfaceDerivsAlgAuv(n,p,U,m,q,V,P,W,xi,eta,d)
    Denom = wders[0][0]
    Denom_du = wders[1][0]
    Denom_dv = wders[0][1]
    dR_du = np.zeros((ncpxi,ncpeta))
    dR_dv = np.zeros((ncpxi,ncpeta))
    #---------------Loop over Non vanishing NURBS Basis Functions-------------#
    #***************Forgot to multiply with control points weights*************# Have to do it
    for ii, ii_value in enumerate(NVu):
        for jj, jj_value in enumerate(NVv):
            BFu = BasisFuns(uspan,xi,p,U)
            BFv = BasisFuns(vspan,eta,q,V)
            # Num = BFu[ii]*BFv[jj]
            Num = BFu[ii]*BFv[jj]**W[ii][jj]
            DBFu = DersBasisFuns(uspan,xi,p,d,U)
            DBFv = DersBasisFuns(vspan,eta,q,d,V)
            Num_du = DBFu[1][ii]*BFv[jj]*W[ii][jj]
            Num_dv = BFu[ii]*DBFv[1][jj]*W[ii][jj]
            #Num_du = DBFu[1][ii]*BFv[jj]
            #Num_dv = BFu[ii]*DBFv[1][jj]
            dR_du[ii_value][jj_value] = Num_du/Denom - Denom_du*Num/(Denom*Denom)
            dR_dv[ii_value][jj_value] = Num_dv/Denom - Denom_dv*Num/(Denom*Denom)
    #---------Flatten (Convert 2D to 1D array) DervsNURBS Function------------#
    #dR_du(0,0)....dR_du(0,1),dR_du(1,0),.....dR_du(1,4).....dR_du(4,0),..........dR_du(4,4)
    #print(dR_du)
    #fdR_du = dR_du.flatten()
    fdR_du = (np.transpose(dR_du)).flatten() #************Cross check this******************#
    #print(fdR_du)
    #fdR_dv = dR_dv.flatten()
    fdR_dv = (np.transpose(dR_dv)).flatten()
    #print(fdR_dv)
    #-------------------------Bu Matrix --------------------------#
    for i2 in range(necp):
        j1= 2*i2
        j2= 2*i2+1
        Bu[0,j1] = fdR_du[i2]
        #print(j1,j2)
        #print(Bu[0,j1])
        Bu[1,j2] = fdR_dv[i2]
        #print(Bu[1,j2])
        Bu[2,j1] = fdR_dv[i2]
        #print(Bu[2,j1])
        Bu[2,j2] = fdR_du[i2]
        #print(Bu[2,j2])

        Be[0,i2] = fdR_du[i2]
        Be[1,i2] = fdR_dv[i2]
    return Bu,Be




def elementRoutine(U_e, T_m):
    Bu=np.zeros((3,nudof*necp))
    Be=np.zeros((2,necp))

    Kt_e = np.zeros(((nudof+nedof)*necp,(nudof+nedof)*necp))

    K_MM=np.zeros((nudof*necp,nudof*necp))
    K_ME=np.zeros((nudof*necp,nedof*necp))
    K_EM=np.zeros((nedof*necp,nudof*necp))
    K_EE=np.zeros((nedof*necp,nedof*necp))
    Fu_int_e=np.zeros((nudof*necp,1))
    Fe_int_e=np.zeros((nedof*necp,1))
    F_int_e=np.zeros(((nudof+nedof)*necp,1))
    sigma_ig=np.zeros((np.shape(GPs_Ws)[0],3,1))
    epsilon_ig=np.zeros((np.shape(GPs_Ws)[0],3,1))
    #F_ext_e=np.zeros(((nudof+nedof)*necp,1)

    #-----------Looping over gauss point-----------#

    for j in range(np.shape(GPs_Ws)[0]):
        file=open('Random.txt','w')
        file.write('Hello')

                #$$$$$$$Have to decide where to place it$$$$$$$$$#
        elU = np.array([0,1]) #Manually defined have to generate using connectivity functions
        elV = np.array([0,1]) #Manually defined have to generate using connectivity functions

        gp = GPs_Ws[j,0:2]
        wg = GPs_Ws[j,2]
        ximas = gp[0]  #Gauss points in Master space
        etamas = gp[1]
        xi = 0.5*((elU[1]-elU[0])*ximas + (elU[1]+elU[0]))     #Paramteric co-ordinates
        eta = 0.5*((elV[1]-elV[0])*etamas + (elV[1]+elV[0]))

        #------Calculating J1 and J2 determinants-----------#

        J1det,J2det = Jacobian12(j,xi,eta,elU,elV)

        #-------------------------Bu Matrix --------------------------#

        Bumatrix,Bematrix = B_matrix(xi,eta,Bu,Be)
        #print('Bu:',Bumatrix)
        #print('Be:',Bematrix)

        # U_u contains mechanical displcement from node 1 to 4, and U_phi contains electric potential from node 1 to 4
        U_u=U_e[0:8]
        U_phi = U_e[8:]

        epsilon = np.matmul(Bumatrix,U_u)
        electric_field = -np.matmul(Bematrix,U_phi)
        epsilon_ig[j]=epsilon
        #print(epsilon)
        #------------------------- C Matrix--------------------------#
        #C=np.array([[2*MU+lamda,lamda,0],[lamda,2*MU+lamda,0],[0,0,MU]])
        #sigma = np.matmul(C,epsilon)
        C, e, k, sigma, Electrical_Displacement = materialRoutine(epsilon,electric_field, T_m)
        sigma_ig[j] = sigma
        #print(C)
        #C=np.array([[139000,74280,0],[74280,115400,0],[0,0,115400]])
        #-------------------------Local Stiffness matrix Ke-------------------#
        CBu=np.matmul(C,Bumatrix)
        BuCBu = np.matmul(np.transpose(Bumatrix),CBu)

        eBe=np.matmul(np.transpose(e),Bematrix)
        BueBe=np.matmul(np.transpose(Bumatrix),eBe)

        eBu=np.matmul(e,Bumatrix)
        BeeBu = np.matmul(np.transpose(Bematrix),eBu)

        kBe=np.matmul(k,Bematrix)
        BekBe = -np.matmul(np.transpose(Bematrix),kBe)

        #------------------Numerical Integration-----------------------#

        #$$$$ Have to add thickness of the plate $$$$$#
        #print('KEE:',K_EE)

        K_MM = K_MM + BuCBu*J1det*J2det*wg*Thick
        K_ME = K_ME + BueBe*J1det*J2det*wg*Thick
        K_EM = K_EM + BeeBu*J1det*J2det*wg*Thick
        K_EE = K_EE + BekBe*J1det*J2det*wg*Thick
        #print('KEE:',K_EE)

        #Arranging to Kt_e matrix 
        Kt_e[0:8,0:8]   = K_MM
        Kt_e[0:8,8:12]  = K_ME
        Kt_e[8:12,0:8]  = K_EM
        Kt_e[8:12,8:12] = K_EE
        #print('K_MM',Kt_e[0:8,0:8])
        #print('K_ME',Kt_e[0:8,8:12])
        #print('K_EM',Kt_e[8:12,0:8])
        #print('K_EE',Kt_e[8:12,8:12])

        #?????? Is this correct way of defining ??????#
        #print('Fu_int_e',Fu_int_e)
        Fu_int_e= Fu_int_e+np.matmul(np.transpose(Bumatrix),sigma)*J1det*J2det*wg*Thick
        Fe_int_e= Fe_int_e+np.matmul(np.transpose(Bematrix),Electrical_Displacement)*J1det*J2det*wg*Thick
        #print('Fu_int_e',Fu_int_e)
        #print(Fe_int_e)

        #Arranging to F_int matrix 
        F_int_e[0:nudof*necp] = Fu_int_e                              # 0 1 2 3 4 5 6 7
        F_int_e[nudof*necp:(nudof*necp+nedof*necp)] = Fe_int_e        # 8 9 10 11
        #print('Fu_int_e',Fu_int_e)
        #print('Fe_int_e',Fe_int_e)
        #print('F_int_e',F_int_e)
        F_ext_e = np.zeros_like(F_int_e)    

    #$$$$$     Didnt return gauss points co-ordinates     $$$$$#    
    #return Kt_e, F_int_e, F_ext_e,sigma
    #print('Stiffness_Matrix:',Kt_e)
    #print('Det of Stiffness_Matrix:',np.linalg.det(Kt_e))
    return Kt_e, F_int_e, F_ext_e, sigma_ig, Electrical_Displacement,epsilon_ig


#----------------Test case------------------#

# U_e = np.zeros((necp*nudof,1))
# Kt_e, F_int_e, F_ext_e,sigma = elementRoutine(U_e, 1)

# print(Kt_e)
# print(F_int_e)
# print(F_ext_e)
# print(sigma)